
import json
import re
import os
import asyncio
from langchain_core.messages import AIMessage, ToolMessage
from langchain_core.prompts import ChatPromptTemplate
from .models import PaperMetadata, ResearchState
from .tools import tools, vectorstore_ref
from langchain_groq import ChatGroq

llm = ChatGroq(
    temperature=0,
    model_name="llama3-70b-8192",
    groq_api_key=os.environ["GROQ_API_KEY"],
    max_retries=3,
    request_timeout=60,
)

agent_prompt_template = ChatPromptTemplate.from_messages([
    ("system", """You are a specialized biomedical research assistant. 
Your goal is to conduct a literature review on a given topic.

You will proceed through the following stages:
1.  **search**: Use the `pubmed_search` and `arxiv_search` tools to find 
relevant papers...
2.  **process**: For promising papers found in the search that have a 
direct PDF URL (e.g., from arXiv), use the `download_and_extract_pdf` 
tool...
3.  **synthesize**: Once you have gathered sufficient abstracts and 
processed relevant full texts..., you will stop using tools. You will 
then analyze the collected abstracts and full texts to write a 
comprehensive summary and identify research gaps.

Instructions:
- Base your actions *only* on the provided `messages` history.
- The user's query is in the first human message.
- Your last message should be the final answer, not a tool call.
- The current stage of research is: **{current_stage}**.
- After using search tools, review the results. If there are relevant 
arXiv papers, consider calling `download_and_extract_pdf` for them.
- Once you believe you have enough information or cannot get more 
full-text, provide your final synthesis without calling any more tools."""),
    ("placeholder", "{messages}"),
])


def agent_node(state: ResearchState) -> dict:
    agent_chain = agent_prompt_template | llm.bind_tools(tools)
    result_message = agent_chain.invoke(state)
    new_search_results = state.get("search_results", [])[:]
    new_processed_papers = state.get("processed_papers", [])[:]
    new_current_stage = state.get("current_stage", "search")

    for msg in state["messages"]:
        if isinstance(msg, ToolMessage):
            try:
                tool_output_data = json.loads(msg.content)
                if isinstance(tool_output_data, list):
                    for item_data in tool_output_data:
                        if isinstance(item_data, dict):
                            try:
                                paper = PaperMetadata(**item_data)
                                if paper.source != "Error":
                                    if paper.paper_id not in [
                                        p.paper_id for p in new_search_results
                                    ]:
                                        new_search_results.append(paper)
                            except Exception:
                                continue
                elif isinstance(tool_output_data, dict):
                    try:
                        paper = PaperMetadata(**tool_output_data)
                        if paper.has_full_text and paper.source != "Error":
                            if paper.paper_id not in [
                                p.paper_id for p in new_processed_papers
                            ]:
                                new_processed_papers.append(paper)
                    except Exception:
                        continue
            except Exception:
                continue

    if getattr(result_message, 'tool_calls', None):
        if any(tc.get('function', {}).get('name') in [
            "pubmed_search", "arxiv_search"
        ] for tc in result_message.tool_calls):
            new_current_stage = "search"
        elif any(tc.get('function', {}).get('name') == 
                 "download_and_extract_pdf" for tc in result_message.tool_calls):
            new_current_stage = "process"

    return {
        "messages": [result_message],
        "search_results": new_search_results,
        "processed_papers": new_processed_papers,
        "current_stage": new_current_stage
    }


def synthesize_node(state: ResearchState) -> dict:
    retriever = vectorstore_ref.as_retriever(search_kwargs={"k": 10})
    retrieved_docs = retriever.invoke(state['research_query'])
    context_text = "\n\n---\n\n".join(
        f"Source URL: {doc.metadata.get('url', 'N/A')}\n\n"
        f"Content: {doc.page_content}"
        for doc in retrieved_docs
    )
    processed_abstracts = ""
    retrieved_urls = {doc.metadata.get('url') for doc in retrieved_docs}
    for paper in state["processed_papers"]:
        if paper.url not in retrieved_urls:
            processed_abstracts += (
                f"\n\nTitle: {paper.title}\nAbstract: {paper.abstract}\n"
                f"URL: {paper.url}"
            )
    if processed_abstracts:
        context_text += (
            f"\n\n--- Processed Paper Abstracts (not in full-text RAG) ---\n"
            f"{processed_abstracts}"
        )

    synthesis_prompt = ChatPromptTemplate.from_template("""
Based on the following research context, please provide a comprehensive 
summary of the findings regarding the user's query: "{query}".

Also, identify and list key research gaps that emerge from the literature.

**Context from Literature:**
{context}

**Output Format:**

### Comprehensive Summary
[Your detailed synthesis of the findings based *only* on the provided 
context.]

### Identified Research Gaps
1. [Clearly stated research gap 1]
2. [Clearly stated research gap 2]
...
""")
    synthesis_chain = synthesis_prompt | llm
    result = synthesis_chain.invoke({
        "query": state['research_query'],
        "context": context_text if context_text else (
            "No relevant full-text documents were processed, and no relevant "
            "abstracts were provided from processed papers."
        )
    })
    summary_content = result.content
    gaps = []
    for line in summary_content.splitlines():
        if re.match(r"^\d+\.\s.*", line.strip()):
            gaps.append(line.strip())
    return {
        "synthesized_summary": summary_content,
        "identified_gaps": gaps,
        "current_stage": "complete",
        "messages": state["messages"] + [AIMessage(content=summary_content)]
    }

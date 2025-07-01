
import json
import re
import os
import asyncio
import uuid
from langchain_core.messages import (
    AIMessage, ToolMessage, HumanMessage
)
from langchain_core.prompts import ChatPromptTemplate
from .models import PaperMetadata, ResearchState
from .tools import tools, vectorstore_ref
from langchain_groq import ChatGroq

llm = ChatGroq(
    temperature=0,
    model_name="llama3-8b-8192",
    groq_api_key=os.environ["GROQ_API_KEY"],
    max_retries=3,
    request_timeout=60,
)

agent_prompt_template = ChatPromptTemplate.from_messages([
    ("system", """You are a specialized biomedical research assistant. 
Your goal is to conduct a literature review on a given topic.

You will proceed through the following stages:
1.  **search**: Use the `pubmed_search`, `arxiv_search`, and `tavily_search` 
tools to find relevant papers and web information.
2.  **process**: For promising papers found in the search, use the 
`download_and_extract_pdf` tool to get their full text content. You should 
only do this if the paper's URL points to a PDF.
3.  **synthesize**: Once papers have been processed, you will stop using 
tools. You will analyze the collected abstracts and full texts to write a 
comprehensive summary and identify research gaps.

Instructions:
- Base your actions *only* on the provided `messages` history.
- The user's query is in the first human message.
- Your last message should be the final answer, not a tool call.
- The current stage of research is: **{current_stage}**."""),
    ("placeholder", "{messages}"),
])


async def agent_node(state: ResearchState):
    """
    The primary node that drives the research process. It invokes the LLM
    to decide the next action based on the current state.
    Handles tool calling explicitly.
    """
    if state["current_stage"] == "search":
        # Generate tool calls for initial search
        tool_calls = []
        
        # Always include PubMed and arXiv
        tool_calls.extend([
            {
                "id": str(uuid.uuid4()),
                "name": "pubmed_search",
                "args": {"query": state["research_query"]}
            },
            {
                "id": str(uuid.uuid4()),
                "name": "arxiv_search",
                "args": {"query": state["research_query"]}
            }
        ])
        
        # Add Tavily if available
        if any(tool.__name__ == "tavily_search" for tool in tools):
            tool_calls.append({
                "id": str(uuid.uuid4()),
                "name": "tavily_search",
                "args": {"query": state["research_query"]}
            })
        
        response = AIMessage(content="", tool_calls=tool_calls)
        
    elif state["current_stage"] == "process":
        # Find arXiv papers that haven't been processed
        papers_to_process = [
            p for p in state["search_results"]
            if (p.source == "arXiv" and p.url and 
                p.paper_id not in [
                    pp.paper_id for pp in state["processed_papers"]
                ])
        ]
        
        if papers_to_process:
            tool_calls = []
            for paper in papers_to_process:
                tool_calls.append({
                    "id": str(uuid.uuid4()),
                    "name": "download_and_extract_pdf",
                    "args": {"url": paper.url, "paper_id": paper.paper_id}
                })
            response = AIMessage(content="", tool_calls=tool_calls)
        else:
            response = AIMessage(
                content="No more arXiv papers to process. Moving to synthesis."
            )
    else:
        # For synthesize or complete stage, use the LLM normally
        agent_chain = agent_prompt_template | llm
        response = await agent_chain.ainvoke(state)

    return {"messages": [response], "current_stage": state["current_stage"]}


async def process_node(state: ResearchState) -> dict:
    """
    A custom node to process tool outputs and update the state.
    This node runs after the 'tools' node.
    """
    new_search_results = state.get("search_results", [])[:]
    new_processed_papers = state.get("processed_papers", [])[:]
    
    # Process recent tool messages
    for msg in reversed(state["messages"]):
        if isinstance(msg, ToolMessage):
            try:
                tool_output_data = json.loads(msg.content)
                if isinstance(tool_output_data, list):
                    for item_data in tool_output_data:
                        if isinstance(item_data, dict):
                            try:
                                # Handle regular paper metadata
                                if "title" in item_data and "abstract" in item_data:
                                    paper = PaperMetadata(**item_data)
                                    if (paper.source != "Error" and 
                                        paper.paper_id not in [
                                            p.paper_id for p in new_search_results
                                        ]):
                                        new_search_results.append(paper)
                                # Handle Tavily results
                                elif ("title" in item_data and "content" in item_data 
                                      and "url" in item_data):
                                    paper = PaperMetadata(
                                        title=item_data.get("title", "No Title"),
                                        abstract=item_data.get("content", "No Content"),
                                        url=item_data.get("url", ""),
                                        year=None,
                                        source="Tavily",
                                        paper_id=item_data.get("url", ""),
                                        has_full_text=False,
                                        authors=[]
                                    )
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
            except json.JSONDecodeError:
                continue
            except Exception:
                continue

    # Determine next stage
    papers_to_process = [
        p for p in new_search_results
        if (p.source == "arXiv" and p.url and 
            p.paper_id not in [pp.paper_id for pp in new_processed_papers])
    ]
    
    next_stage = "process" if papers_to_process else "synthesize"

    return {
        "search_results": new_search_results,
        "processed_papers": new_processed_papers,
        "current_stage": next_stage,
        "messages": state["messages"]
    }


def synthesize_node(state: ResearchState) -> dict:
    """
    Generates the final summary and identifies research gaps based on all
    collected data. Uses both abstracts and processed full texts.
    """
    # Get context from processed PDFs
    retriever = vectorstore_ref.as_retriever(search_kwargs={"k": 10})
    retrieved_docs = retriever.invoke(state['research_query'])
    
    combined_context = ""
    
    if retrieved_docs:
        combined_context += "**Full Text Content (from Processed PDFs):**\n\n"
        combined_context += "\n\n---\n\n".join(
            f"Source URL: {doc.metadata.get('url', 'N/A')}\n\n"
            f"Content: {doc.page_content[:1000]}..."
            for doc in retrieved_docs
        )
        combined_context += "\n\n"

    # Add abstracts from search results
    processed_paper_ids = {p.paper_id for p in state["processed_papers"]}
    search_results_to_include = [
        p for p in state["search_results"]
        if (p.paper_id not in processed_paper_ids and p.abstract and 
            p.abstract not in ["No Abstract Found.", "No Content"])
    ]

    if search_results_to_include:
        combined_context += "**Context from Initial Search (Abstracts/Snippets):**\n\n"
        combined_context += "\n\n---\n\n".join(
            f"Source: {p.source}, Title: {p.title}\n\nContent: {p.abstract}"
            for p in search_results_to_include
        )
        combined_context += "\n\n"

    synthesis_prompt = ChatPromptTemplate.from_template("""
Based on the following research context, which includes both full text content 
(where available) and abstracts/snippets from initial searches, please provide 
a comprehensive summary of the findings regarding the user's query: "{query}".

Prioritize information from the full text content when available.

Also, identify and list key research gaps that emerge from the literature.

**Research Context:**
{context}

**Output Format:**

### Comprehensive Summary
[Your detailed synthesis of the findings based *only* on the provided context.]

### Identified Research Gaps
1. [Clearly stated research gap 1]
2. [Clearly stated research gap 2]
...
""")
    
    synthesis_chain = synthesis_prompt | llm
    result = synthesis_chain.invoke({
        "query": state['research_query'],
        "context": combined_context if combined_context else (
            "No relevant documents or abstracts were found or processed."
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

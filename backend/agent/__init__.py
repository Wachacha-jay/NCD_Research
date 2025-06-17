
from langchain_core.messages import HumanMessage
from .env_setup import *  # noqa: F403,F401
from .models import ResearchState  # noqa: F401
from .workflow import build_workflow

# Provide direct reference to tools/vectorstore in tools.py
from .tools import vectorstore_ref  # noqa: F401

_agent = None


async def run_research_agent(query: str):
    global _agent
    if _agent is None:
        _agent = build_workflow()
    initial_state = {
        "messages": [HumanMessage(content=query)],
        "research_query": query,
        "search_results": [],
        "processed_papers": [],
        "synthesized_summary": "",
        "identified_gaps": [],
        "current_stage": "search",
        "error": None
    }
    final_output = {}
    try:
        final_state = await _agent.ainvoke(
            initial_state, {"recursion_limit": 100}
        )
        final_output = {
            "query": final_state["research_query"],
            "summary": final_state["synthesized_summary"],
            "gaps": final_state["identified_gaps"],
            "search_results_count": len(final_state["search_results"]),
            "processed_papers_count": len(final_state["processed_papers"]),
            "status": "completed",
            "error": final_state["error"]
        }
    except Exception as e:
        final_output = {
            "query": query,
            "summary": "An error occurred during research.",
            "gaps": [],
            "search_results_count": 0,
            "processed_papers_count": 0,
            "status": "failed",
            "error": str(e)
        }
    return final_output

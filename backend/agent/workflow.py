
from langgraph.graph import StateGraph, END
from langgraph.prebuilt import ToolNode
from typing import Literal
from .models import ResearchState
from .tools import tools
from .nodes import agent_node, synthesize_node


def should_continue(state: ResearchState) -> Literal["tools", "synthesize", "end"]:
    last_message = state["messages"][-1]
    from langchain_core.messages import AIMessage
    if isinstance(last_message, AIMessage) and not getattr(
        last_message, 'tool_calls', None
    ):
        return "synthesize"
    return "tools"


def build_workflow():
    workflow = StateGraph(ResearchState)
    tool_node = ToolNode(tools)
    workflow.add_node("agent", agent_node)
    workflow.add_node("tools", tool_node)
    workflow.add_node("synthesize", synthesize_node)
    workflow.set_entry_point("agent")
    workflow.add_conditional_edges(
        "agent",
        should_continue,
        {
            "tools": "tools",
            "synthesize": "synthesize",
            "end": END
        }
    )
    workflow.add_edge("tools", "agent")
    workflow.add_edge("synthesize", END)
    return workflow.compile()

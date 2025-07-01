
from langgraph.graph import StateGraph, END
from langgraph.prebuilt import ToolNode
from typing import Literal
from .models import ResearchState
from .tools import tools
from .nodes import agent_node, synthesize_node, process_node


def should_continue(state: ResearchState) -> Literal["tools", "process", "synthesize", "end"]:
    """Conditional logic to route the workflow."""
    last_message = state["messages"][-1]
    
    # If the last message is a tool call, go to tools
    if hasattr(last_message, 'tool_calls') and last_message.tool_calls:
        return "tools"
    
    # If the current stage is search and no tool calls, move to process
    if state["current_stage"] == "search":
        return "process"
    
    # If the current stage is synthesize, end
    if state["current_stage"] == "synthesize":
        return "end"
    
    # Default to end
    return "end"


def build_workflow():
    workflow = StateGraph(ResearchState)
    tool_node = ToolNode(tools)
    
    # Add nodes
    workflow.add_node("agent", agent_node)
    workflow.add_node("tools", tool_node)
    workflow.add_node("process", process_node)
    workflow.add_node("synthesize", synthesize_node)
    
    # Set entry point
    workflow.set_entry_point("agent")
    
    # Add conditional edges from agent
    workflow.add_conditional_edges(
        "agent",
        should_continue,
        {
            "tools": "tools",
            "process": "process",
            "synthesize": "synthesize",
            "end": END
        }
    )
    
    # From tools, always go to process
    workflow.add_edge("tools", "process")
    
    # From process, decide whether to go back to agent or synthesize
    workflow.add_conditional_edges(
        "process",
        lambda state: "agent" if state["current_stage"] == "process" else "synthesize",
        {"agent": "agent", "synthesize": "synthesize"}
    )
    
    # From synthesize, always end
    workflow.add_edge("synthesize", END)
    
    return workflow.compile()


from typing import List, Optional, Literal, Annotated, TypedDict
from langchain_core.messages import BaseMessage
from pydantic import BaseModel, Field

class PaperMetadata(BaseModel):
    title: str
    abstract: str
    url: str
    year: Optional[str] = None
    source: Literal["PubMed", "arXiv", "Semantic Scholar", "PDF", "Error"]
    paper_id: str
    has_full_text: bool = False
    authors: List[str] = Field(default_factory=list)

class ResearchState(TypedDict):
    messages: Annotated[List[BaseMessage], lambda x, y: x + y]
    research_query: str
    search_results: Annotated[List[PaperMetadata], lambda x, y: x + y]
    processed_papers: Annotated[List[PaperMetadata], lambda x, y: x + y]
    synthesized_summary: str
    identified_gaps: Annotated[List[str], lambda x, y: list(set(x + y))]
    current_stage: Literal["search", "process", "synthesize", "complete"]
    error: Optional[str]

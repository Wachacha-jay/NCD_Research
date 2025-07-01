
import asyncio
import json
import uuid
from langchain_core.tools import tool
from .models import PaperMetadata
from Bio import Entrez
import arxiv
import aiohttp
from PyPDF2 import PdfReader
from langchain_core.documents import Document
from io import BytesIO
import os

from langchain_community.embeddings import HuggingFaceEmbeddings
from langchain_community.vectorstores import Chroma

# Try to import Tavily
try:
    from langchain_tavily import TavilySearch
    TAVILY_AVAILABLE = True
except ImportError:
    TAVILY_AVAILABLE = False
    print("Tavily not available - install with: pip install langchain-tavily")

CHROMA_DB_PATH = "./chroma_db_agent"
os.makedirs(CHROMA_DB_PATH, exist_ok=True)
embeddings = HuggingFaceEmbeddings(model_name="BAAI/bge-small-en-v1.5")
vectorstore = Chroma(
    embedding_function=embeddings, persist_directory=CHROMA_DB_PATH
)


@tool
async def pubmed_search(query: str, limit: int = 5):
    """
    Searches PubMed for biomedical literature based on a query.
    Returns a list of paper metadata.
    """
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=limit)
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        papers = []
        if id_list:
            fetch_handle = Entrez.efetch(
                db="pubmed", id=",".join(id_list), retmode="xml"
            )
            data = Entrez.read(fetch_handle)
            fetch_handle.close()
            for article in data.get("PubmedArticle", []):
                medline = article.get("MedlineCitation", {})
                article_info = medline.get("Article", {})
                if not article_info:
                    continue
                abstract_parts = article_info.get(
                    "Abstract", {}
                ).get("AbstractText", [])
                paper = PaperMetadata(
                    title=str(article_info.get(
                        "ArticleTitle", "No Title Found"
                    )),
                    abstract=" ".join(abstract_parts) if abstract_parts else (
                        "No Abstract Found."
                    ),
                    url=f"https://pubmed.ncbi.nlm.nih.gov/"
                        f"{medline.get('PMID', '')}/",
                    year=article_info.get("Journal", {}).get(
                        "JournalIssue", {}
                    ).get("PubDate", {}).get("Year"),
                    source="PubMed",
                    paper_id=str(medline.get("PMID", "")),
                    authors=[
                        f"{a.get('LastName', '')}, {a.get('ForeName', '')}"
                        for a in article_info.get("AuthorList", [])
                    ]
                )
                papers.append(paper)
        return papers
    except Exception as e:
        return [PaperMetadata(
            title="PubMed Search Error",
            abstract=f"An error occurred: {str(e)}",
            url="",
            source="Error",
            paper_id="error-pubmed"
        )]


@tool
async def arxiv_search(query: str, limit: int = 5):
    """
    Searches arXiv for pre-print papers, especially in physics, math, and
    computer science. Returns a list of paper metadata.
    """
    try:
        client = arxiv.Client()
        search = arxiv.Search(
            query=query,
            max_results=limit,
            sort_by=arxiv.SortCriterion.Relevance
        )
        papers = []
        for result in await asyncio.to_thread(
            lambda: list(client.results(search))
        ):
            paper = PaperMetadata(
                title=result.title,
                abstract=result.summary.replace("\n", " "),
                url=result.pdf_url,
                year=str(result.published.year),
                source="arXiv",
                paper_id=result.entry_id.split('/')[-1],
                authors=[author.name for author in result.authors],
                has_full_text=True
            )
            papers.append(paper)
        return papers
    except Exception as e:
        return [PaperMetadata(
            title="arXiv Search Error",
            abstract=f"An error occurred: {str(e)}",
            url="",
            source="Error",
            paper_id="error-arxiv"
        )]


@tool
async def download_and_extract_pdf(url: str, paper_id: str):
    """
    Downloads a PDF from a URL, extracts its text, and stores it in a
    vector database. This tool should be used when a paper's full text
    is needed for synthesis.
    """
    headers = {"User-Agent": "ResearchAgent/1.0"}
    async with aiohttp.ClientSession() as session:
        try:
            async with session.get(url, headers=headers, timeout=30) as resp:
                resp.raise_for_status()
                content = await resp.read()
                reader = PdfReader(BytesIO(content))
                text = "\n".join(
                    page.extract_text() or "" for page in reader.pages
                )
                if not text:
                    raise ValueError(
                        "PDF text extraction resulted in empty content."
                    )
                doc = Document(
                    page_content=text,
                    metadata={
                        "paper_id": paper_id, "source": "PDF", "url": url
                    }
                )
                await vectorstore.aadd_documents([doc])
                return PaperMetadata(
                    title="Successfully Extracted PDF",
                    abstract=text[:500] + "..." if len(text) > 500 else text,
                    url=url,
                    paper_id=paper_id,
                    has_full_text=True,
                    source="PDF"
                )
        except Exception as e:
            return PaperMetadata(
                title="PDF Processing Error",
                abstract=str(e),
                url=url,
                paper_id=paper_id,
                has_full_text=False,
                source="Error"
            )


# Initialize tools list
tools = [pubmed_search, arxiv_search, download_and_extract_pdf]

# Add Tavily search if available
if TAVILY_AVAILABLE:
    tavily_search = TavilySearch(
        max_results=5,
        tavily_api_key=os.getenv("TAVILY_API_KEY")
    )
    tools.append(tavily_search)

vectorstore_ref = vectorstore

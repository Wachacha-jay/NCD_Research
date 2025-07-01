
import os
from Bio import Entrez

# Environment variable config
if "GROQ_API_KEY" not in os.environ:
    os.environ["GROQ_API_KEY"] = "YOUR_GROQ_API_KEY"
if "LANGCHAIN_API_KEY" not in os.environ:
    os.environ["LANGCHAIN_API_KEY"] = "YOUR_LANGCHAIN_API_KEY"
    os.environ["LANGCHAIN_TRACING_V2"] = "true"
if "PUBMED_API_KEY" not in os.environ:
    os.environ["PUBMED_API_KEY"] = ""
if "TAVILY_API_KEY" not in os.environ:
    os.environ["TAVILY_API_KEY"] = "YOUR_TAVILY_API_KEY"

# Best practice: provide an email for NCBI
Entrez.email = os.getenv("ENTREZ_EMAIL", "your.email@example.com")
if os.getenv("PUBMED_API_KEY"):
    Entrez.api_key = os.getenv("PUBMED_API_KEY")

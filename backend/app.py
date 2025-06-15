
import os
import asyncio
from flask import Flask, request, jsonify
from flask_cors import CORS

# Import your full agent logic here.
# For simplicity, we'll assume your research agent code (classes, agent_node, run_research_agent, etc.)
# is pasted into this file above this API section, or better: moved into a backend/agent.py module.

# Let's assume all agent code is above/below or imported:
#
# from agent import run_research_agent

# ------------------------
# FLASK APP/API ENDPOINTS
# ------------------------

app = Flask(__name__)
CORS(app)  # (optional) allow requests from React frontend

# -- ASYNC EVENT LOOP SUPPORT FOR FLASK --
def run_async(func):
    def wrapper(*args, **kwargs):
        return asyncio.run(func(*args, **kwargs))
    wrapper.__name__ = func.__name__
    return wrapper

@app.route('/api/research', methods=['POST'])
def research_handler():
    data = request.get_json()
    query = data.get("query")
    if not query:
        return jsonify({"status": "error", "message": "Missing 'query'"}), 400
    try:
        # Run the async agent with a new event loop
        result = asyncio.run(run_research_agent(query))
        if result["status"] == "completed":
            return jsonify({
                "status": "success",
                "summary": result["summary"],
                "gaps": result["gaps"],
                "message": "Literature analysis completed successfully.",
                "query": result["query"]
            })
        else:
            return jsonify({
                "status": "failed",
                "message": result.get("error", "Unknown error."),
                "summary": "",
                "gaps": [],
                "query": result["query"]
            }), 500
    except Exception as e:
        return jsonify({"status": "error", "message": str(e)}), 500

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000)

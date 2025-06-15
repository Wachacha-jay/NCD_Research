
import os
import asyncio
from flask import Flask, request, jsonify
from flask_cors import CORS

# Import the run_research_agent function from the new modular agent
from agent import run_research_agent

app = Flask(__name__)
CORS(app)

@app.route('/api/research', methods=['POST'])
def research_handler():
    data = request.get_json()
    query = data.get("query")
    if not query:
        return jsonify({"status": "error", "message": "Missing 'query'"}), 400
    try:
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

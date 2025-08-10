import os
import requests

# Gemini API key from environment
API_KEY = "AIzaSyAXDfE1yr7viaqwrgVY_JGFnrmD_FQB_Vo"
MODEL = "gemini-1.5-flash"  # You can change to gemini-1.5-pro for more capability

SYSTEM_PROMPT = """You are BioPharm-Guide, a concise, factual assistant for
questions about proteins, amino acids, enzyme function, pathways, structural biology,
binding/affinity concepts, pharmacology, PK/PD, and general drug mechanisms.
Stay rigorous and cite standard knowledge (e.g., “common textbook consensus”) when unsure.
Avoid medical advice; provide informational content only. Keep answers compact unless asked for depth.
If a question is outside this scope, say so briefly and suggest a related bio topic. Answer in 2-3 sentences, followed by an empty line, following by your brief citation in this format: Source: <brief source info>."""

# Initialize persistent chat history
history = [{"role": "system", "content": SYSTEM_PROMPT}]

def ChatBot(question: str, context: str | None = None) -> str:
    """Send a question to the BioPharm-Guide chatbot and return its answer as a string.
    Optionally include 'context' that will be prepended to the user message.
    """
    if context:
        message = f"Context for reference (use if relevant):\n{context}\n\nQuestion: {question}"
    else:
        message = question
    history.append({"role": "user", "content": message})

    # Prepare messages for Gemini API
    contents = []
    for msg in history:
        role = msg["role"]
        if role == "system":
            contents.append({"role": "user", "parts": [{"text": msg["content"]}]})
        elif role == "user":
            contents.append({"role": "user", "parts": [{"text": msg["content"]}]})
        elif role == "assistant":
            contents.append({"role": "model", "parts": [{"text": msg["content"]}]})

    # If API key isn't available, provide a deterministic, offline fallback answer
    if not API_KEY:
        fallback = "This is an offline informational reply. " \
                   "I can summarize provided context and general protein concepts, " \
                   "but I’m not using a live LLM.\n\n"
        if context:
            fallback += "Context summary: Top similarity hits and properties received.\n"
        fallback += "Source: local rules + provided context"
        history.append({"role": "assistant", "content": fallback})
        return fallback

    url = f"https://generativelanguage.googleapis.com/v1beta/models/{MODEL}:generateContent?key={API_KEY}"
    resp = requests.post(url, json={"contents": contents}, timeout=30)

    if resp.status_code != 200:
        # Graceful fallback
        answer = f"LLM unavailable (HTTP {resp.status_code}). " \
                 f"Providing offline summary.\n\nSource: local rules + provided context"
        history.append({"role": "assistant", "content": answer})
        return answer

    data = resp.json()
    try:
        answer = data["candidates"][0]["content"]["parts"][0]["text"].strip()
    except (KeyError, IndexError):
        answer = "LLM returned an unexpected format. Source: local rules + provided context"

    history.append({"role": "assistant", "content": answer})
    return answer


# Example usage
if __name__ == "__main__":
    q = "What is the role of tryptophan in serotonin synthesis?"
    print(ChatBot(q))

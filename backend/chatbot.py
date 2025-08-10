import os
import requests

# Gemini API key
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

def ChatBot(question: str) -> str:
    """Send a question to the BioPharm-Guide chatbot and return its answer as a string."""
    history.append({"role": "user", "content": question})

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

    url = f"https://generativelanguage.googleapis.com/v1beta/models/{MODEL}:generateContent?key={API_KEY}"
    resp = requests.post(url, json={"contents": contents})

    if resp.status_code != 200:
        raise RuntimeError(f"Gemini API error {resp.status_code}: {resp.text}")

    data = resp.json()
    try:
        answer = data["candidates"][0]["content"]["parts"][0]["text"].strip()
    except (KeyError, IndexError):
        raise RuntimeError(f"Unexpected Gemini response format: {data}")

    history.append({"role": "assistant", "content": answer})
    return answer


# Example usage
if __name__ == "__main__":
    q = "What is the role of tryptophan in serotonin synthesis?"
    print(ChatBot(q))

#!/usr/bin/env python3
"""
dev_server.py — Start the FastAPI backend for local development.

Usage (from project root):
    python dev_server.py

This script adds the backend/ directory to sys.path so all non-prefixed
imports (config, db, models, routers, services) resolve correctly.
"""
import subprocess
import sys
from pathlib import Path

BACKEND_DIR = Path(__file__).parent / "backend"

if __name__ == "__main__":
    env_file = BACKEND_DIR / ".env"
    if not env_file.exists():
        example = BACKEND_DIR / ".env.example"
        if example.exists():
            import shutil
            shutil.copy(example, env_file)
            print(f"[dev_server] Created {env_file} from .env.example — please fill in LLM_API_KEY")
        else:
            print("[dev_server] Warning: .env not found and no .env.example to copy from")

    cmd = [
        sys.executable, "-m", "uvicorn", "main:app",
        "--host", "0.0.0.0",
        "--port", "8000",
        "--reload",
        "--reload-dir", str(BACKEND_DIR),
        "--env-file", str(env_file),
    ]
    subprocess.run(cmd, cwd=str(BACKEND_DIR))

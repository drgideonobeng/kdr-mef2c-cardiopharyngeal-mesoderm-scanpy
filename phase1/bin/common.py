"""
common.py

Small shared helpers used by every script in phase1/bin.
Kept deliberately minimal -- a logger and a fallback seed constant.
No scientific logic lives here.
"""

# Fallback used only if a script doesn't receive --seed explicitly.
DEFAULT_SEED = 42


def log(step: str, msg: str) -> None:
    """Print a numbered log line, e.g. log('03', 'filtered 120 cells')."""
    print(f"[{step}] {msg}", flush=True)

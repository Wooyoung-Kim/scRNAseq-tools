import sys

from ..llm import KNOWN_PROVIDERS, has_api_key, resolve_model, resolve_provider


def register(subparsers) -> None:
    parser = subparsers.add_parser("doctor", help="Check environment and configuration")
    parser.set_defaults(func=run)


def run(_args) -> int:
    provider = resolve_provider()
    model = resolve_model()

    print(f"Python: {sys.version.split()[0]}")
    print(f"LLM provider: {provider}")
    if model:
        print(f"LLM model: {model}")

    if provider == "none":
        print("LLM: disabled (default)")
        return 0

    if provider not in KNOWN_PROVIDERS:
        print("LLM: custom provider (ensure your own API key handling)")
        return 0

    if has_api_key(provider):
        print("LLM: API key detected")
        return 0

    if provider == "openai":
        print("LLM: missing OPENAI_API_KEY")
    elif provider == "anthropic":
        print("LLM: missing ANTHROPIC_API_KEY")
    else:
        print("LLM: provider ready (no API key required)")

    return 0


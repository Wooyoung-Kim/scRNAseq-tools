from .. import __version__
from ..llm import has_api_key, resolve_model, resolve_provider


def register(subparsers) -> None:
    parser = subparsers.add_parser("info", help="Show basic info")
    parser.set_defaults(func=run)


def run(_args) -> int:
    provider = resolve_provider()
    model = resolve_model()
    if provider == "none":
        key_status = "n/a"
    else:
        key_status = "set" if has_api_key(provider) else "missing"

    print(f"scRNAseq-tools {__version__}")
    print(f"LLM provider: {provider}")
    if model:
        print(f"LLM model: {model}")
    print(f"API key: {key_status}")
    return 0


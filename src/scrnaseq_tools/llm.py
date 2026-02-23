import os

KNOWN_PROVIDERS = {"none", "openai", "anthropic", "ollama"}


def resolve_provider(env=os.environ) -> str:
    provider = (env.get("LLM_PROVIDER") or "").strip().lower()
    if provider:
        return provider
    if env.get("OPENAI_API_KEY"):
        return "openai"
    if env.get("ANTHROPIC_API_KEY"):
        return "anthropic"
    return "none"


def resolve_model(env=os.environ) -> str:
    return (
        env.get("LLM_MODEL")
        or env.get("MODEL")
        or env.get("OPENAI_MODEL")
        or env.get("ANTHROPIC_MODEL")
        or ""
    )


def has_api_key(provider: str, env=os.environ) -> bool:
    provider = provider.lower()
    if provider == "openai":
        return bool(env.get("OPENAI_API_KEY"))
    if provider == "anthropic":
        return bool(env.get("ANTHROPIC_API_KEY"))
    if provider == "ollama":
        return True
    return False


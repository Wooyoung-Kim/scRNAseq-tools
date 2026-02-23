import json
import os
import shlex
import subprocess
import sys
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, List, Optional


PROVIDERS = (
    "ollama",
    "vllm",
    "lmstudio",
    "openai",
    "anthropic",
    "codex-cli",
    "claude-cli",
)
DEFAULT_PROVIDER = "claude-cli"


def register(subparsers) -> None:
    parser = subparsers.add_parser("chat", help="Interactive chat with local or remote models")
    parser.add_argument(
        "--provider",
        choices=PROVIDERS,
        default=None,
        help="Provider to use (prompted if omitted)",
    )
    parser.add_argument("--model", default=None, help="Model name")
    parser.add_argument("--base-url", default=None, help="Override API base URL")
    parser.add_argument("--system", default=os.environ.get("SCRNASEQ_SYSTEM_PROMPT", ""), help="System prompt")
    parser.add_argument("--max-tokens", type=int, default=None, help="Max tokens for responses")
    parser.add_argument("--temperature", type=float, default=None, help="Sampling temperature")
    parser.add_argument(
        "--no-history",
        action="store_true",
        help="Do not send conversation history",
    )
    parser.add_argument(
        "--cli-cmd",
        default=None,
        help="Override codex/claude CLI command",
    )
    parser.add_argument(
        "--cli-args",
        default="",
        help="Extra arguments for codex/claude CLI (quoted string)",
    )
    parser.set_defaults(func=run)


def prompt_provider(default: str) -> str:
    options = [
        ("ollama", "Local Ollama"),
        ("vllm", "Local vLLM (OpenAI-compatible)"),
        ("lmstudio", "Local LM Studio (OpenAI-compatible)"),
        ("openai", "OpenAI API key"),
        ("anthropic", "Anthropic API key"),
        ("codex-cli", "codex CLI (login-based)"),
        ("claude-cli", "claude CLI (login-based)"),
    ]
    print("Select provider:")
    for idx, (key, label) in enumerate(options, start=1):
        print(f"  {idx}) {label} [{key}]")
    print(f"Press Enter for default: {default}")
    while True:
        choice = input("> ").strip().lower()
        if not choice:
            return default
        if choice.isdigit():
            index = int(choice) - 1
            if 0 <= index < len(options):
                return options[index][0]
        for key, _label in options:
            if choice == key:
                return key
        print("Invalid selection. Enter a number or provider key.")


def detect_cli_identity(cmd: str) -> Optional[str]:
    if os.environ.get("SCRNASEQ_DISABLE_CLI_DETECT", "").lower() in ("1", "true", "yes"):
        return None
    attempts = [
        [cmd, "--version"],
        [cmd, "-v"],
        [cmd, "--help"],
    ]
    for attempt in attempts:
        try:
            result = subprocess.run(
                attempt,
                capture_output=True,
                text=True,
                timeout=2,
            )
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
        output = (result.stdout + result.stderr).lower()
        if "codex" in output or "openai codex" in output:
            return "codex"
        if "claude" in output or "anthropic" in output:
            return "claude"
    return None


def maybe_bootstrap_agents(project_dir: Path) -> None:
    if os.environ.get("SCRNASEQ_DISABLE_AGENTS_BOOTSTRAP", "").lower() in ("1", "true", "yes"):
        return
    agents_path = project_dir / "AGENTS.md"
    claude_dir = project_dir / ".claude"
    if agents_path.exists() or not claude_dir.is_dir():
        return
    content = "\n".join(
        [
            "# AGENTS",
            "",
            "Use this project's .claude configuration when analyzing.",
            "Read and follow these files before starting work:",
            "- .claude/README.md (project overview and rules)",
            "- .claude/settings.json (tools and permissions)",
            "- .claude/agents/annotation-executor.md (annotation workflow)",
            "",
            "If instructions conflict, prefer the .claude guidance.",
            "",
        ]
    )
    try:
        agents_path.write_text(content, encoding="utf-8")
    except OSError:
        return


def maybe_bootstrap_claude_md(project_dir: Path) -> None:
    if os.environ.get("SCRNASEQ_DISABLE_CLAUDE_BOOTSTRAP", "").lower() in ("1", "true", "yes"):
        return
    claude_md = project_dir / "CLAUDE.md"
    claude_dir = project_dir / ".claude"
    if claude_md.exists() or not claude_dir.is_dir():
        return
    content = "\n".join(
        [
            "# Claude Project Instructions",
            "",
            "This project uses the .claude configuration directory.",
            "Please read and follow:",
            "- .claude/README.md",
            "- .claude/settings.json",
            "- .claude/agents/annotation-executor.md",
            "",
        ]
    )
    try:
        claude_md.write_text(content, encoding="utf-8")
    except OSError:
        return


def post_json(url: str, headers: Dict[str, str], payload: Dict) -> Dict:
    data = json.dumps(payload).encode("utf-8")
    request = urllib.request.Request(url, data=data, headers=headers, method="POST")
    with urllib.request.urlopen(request, timeout=120) as response:
        return json.loads(response.read().decode("utf-8"))


def openai_chat(
    base_url: str,
    model: str,
    messages: List[Dict[str, str]],
    api_key: Optional[str],
    temperature: Optional[float],
    max_tokens: Optional[int],
) -> str:
    url = base_url.rstrip("/") + "/v1/chat/completions"
    headers = {"Content-Type": "application/json"}
    if api_key:
        headers["Authorization"] = f"Bearer {api_key}"
    payload: Dict = {"model": model, "messages": messages}
    if temperature is not None:
        payload["temperature"] = temperature
    if max_tokens is not None:
        payload["max_tokens"] = max_tokens
    response = post_json(url, headers, payload)
    return response["choices"][0]["message"]["content"]


def ollama_chat(
    base_url: str,
    model: str,
    messages: List[Dict[str, str]],
    temperature: Optional[float],
    max_tokens: Optional[int],
) -> str:
    url = base_url.rstrip("/") + "/api/chat"
    payload: Dict = {"model": model, "messages": messages, "stream": False}
    options: Dict[str, float] = {}
    if temperature is not None:
        options["temperature"] = temperature
    if max_tokens is not None:
        options["num_predict"] = max_tokens
    if options:
        payload["options"] = options
    response = post_json(url, {"Content-Type": "application/json"}, payload)
    return response["message"]["content"]


def anthropic_chat(
    base_url: str,
    model: str,
    system_prompt: str,
    messages: List[Dict[str, str]],
    api_key: str,
    temperature: Optional[float],
    max_tokens: Optional[int],
) -> str:
    url = base_url.rstrip("/") + "/v1/messages"
    headers = {
        "Content-Type": "application/json",
        "x-api-key": api_key,
        "anthropic-version": "2023-06-01",
    }
    payload: Dict = {
        "model": model,
        "max_tokens": max_tokens or 1024,
        "messages": messages,
    }
    if system_prompt:
        payload["system"] = system_prompt
    if temperature is not None:
        payload["temperature"] = temperature
    response = post_json(url, headers, payload)
    return response["content"][0]["text"]


def run_cli(cmd: str, extra_args: str) -> int:
    args = shlex.split(extra_args) if extra_args else []
    full_cmd = [cmd] + args
    try:
        return subprocess.call(full_cmd)
    except FileNotFoundError:
        print(f"CLI not found: {cmd}")
        return 1


def run(args) -> int:
    provider = args.provider or prompt_provider(DEFAULT_PROVIDER)

    if provider == "codex-cli":
        cmd = args.cli_cmd or os.environ.get("SCRNASEQ_CODEX_CMD", "codex")
        if not args.cli_cmd:
            identity = detect_cli_identity(cmd)
            if identity == "claude":
                print("Detected Claude CLI behind codex; continuing.")
            elif identity == "codex":
                print("Detected Codex CLI.")
        maybe_bootstrap_agents(Path.cwd())
        return run_cli(cmd, args.cli_args)
    if provider == "claude-cli":
        cmd = args.cli_cmd or os.environ.get("SCRNASEQ_CLAUDE_CMD", "claude")
        maybe_bootstrap_claude_md(Path.cwd())
        if not args.cli_cmd:
            identity = detect_cli_identity(cmd)
            if identity == "codex":
                print("Detected Codex CLI behind claude; bootstrapping AGENTS.md if missing.")
                maybe_bootstrap_agents(Path.cwd())
        return run_cli(cmd, args.cli_args)

    system_prompt = args.system or ""
    history: List[Dict[str, str]] = []
    if system_prompt and provider != "anthropic":
        history.append({"role": "system", "content": system_prompt})

    print("Type /exit to quit, /reset to clear history, /help for commands.")
    while True:
        try:
            user_input = input("you> ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            return 0
        if not user_input:
            continue
        if user_input in ("/exit", "/quit"):
            return 0
        if user_input == "/reset":
            history = []
            if system_prompt and provider != "anthropic":
                history.append({"role": "system", "content": system_prompt})
            print("history cleared")
            continue
        if user_input == "/help":
            print("commands: /exit, /quit, /reset, /help")
            continue

        if args.no_history:
            messages = []
            if system_prompt and provider != "anthropic":
                messages.append({"role": "system", "content": system_prompt})
        else:
            messages = list(history)

        messages.append({"role": "user", "content": user_input})

        try:
            reply = dispatch_chat(provider, args, system_prompt, messages)
        except (urllib.error.HTTPError, urllib.error.URLError, KeyError, ValueError) as exc:
            print(f"error: {exc}")
            continue

        print(f"assistant> {reply}")
        if not args.no_history:
            history.append({"role": "user", "content": user_input})
            history.append({"role": "assistant", "content": reply})


def dispatch_chat(provider: str, args, system_prompt: str, messages: List[Dict[str, str]]) -> str:
    if provider == "ollama":
        base_url = args.base_url or os.environ.get("OLLAMA_HOST", "http://localhost:11434")
        model = args.model or os.environ.get("OLLAMA_MODEL", "llama3.1")
        return ollama_chat(base_url, model, messages, args.temperature, args.max_tokens)

    if provider in ("vllm", "lmstudio", "openai"):
        base_url = args.base_url
        if not base_url:
            if provider == "vllm":
                base_url = os.environ.get("VLLM_BASE_URL", "http://localhost:8000")
            elif provider == "lmstudio":
                base_url = os.environ.get("LMSTUDIO_BASE_URL", "http://localhost:1234")
            else:
                base_url = os.environ.get("OPENAI_BASE_URL", "https://api.openai.com")

        api_key = None
        if provider == "openai":
            api_key = os.environ.get("OPENAI_API_KEY")
            if not api_key:
                raise ValueError("OPENAI_API_KEY is required for OpenAI provider")

        if provider == "openai":
            model = args.model or os.environ.get("OPENAI_MODEL") or os.environ.get("LLM_MODEL") or "gpt-4.1-mini"
        elif provider == "vllm":
            model = args.model or os.environ.get("VLLM_MODEL") or os.environ.get("OPENAI_MODEL") or "local-model"
        else:
            model = args.model or os.environ.get("LMSTUDIO_MODEL") or os.environ.get("OPENAI_MODEL") or "local-model"

        return openai_chat(base_url, model, messages, api_key, args.temperature, args.max_tokens)

    if provider == "anthropic":
        base_url = args.base_url or os.environ.get("ANTHROPIC_BASE_URL", "https://api.anthropic.com")
        api_key = os.environ.get("ANTHROPIC_API_KEY")
        if not api_key:
            raise ValueError("ANTHROPIC_API_KEY is required for Anthropic provider")
        model = args.model or os.environ.get("ANTHROPIC_MODEL") or os.environ.get("LLM_MODEL") or "claude-3-5-sonnet-20240620"
        return anthropic_chat(
            base_url,
            model,
            system_prompt,
            messages,
            api_key,
            args.temperature,
            args.max_tokens,
        )

    raise ValueError(f"Unknown provider: {provider}")

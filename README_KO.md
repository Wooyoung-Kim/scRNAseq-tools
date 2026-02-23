# scRNAseq-tools

단일세포 RNA-seq 워크플로우를 위한 CLI 중심 도구입니다. 기본 모드는 외부 LLM API를 사용하지 않습니다.

## 빠른 시작

```bash
python -m pip install -e .
scRNAseq-tools --help
scRNAseq-tools info
scRNAseq-tools init --target /path/to/project
```

## 단계별 사용 가이드

### 프로젝트 초기화 (.claude만 설치)
1) `python -m pip install -e .`
2) `scRNAseq-tools init --target /path/to/project`
3) `cd /path/to/project` 후 `claude` 또는 `codex` 실행

### 프로젝트 부트스트랩 (MCP 의존성 포함)
1) `scRNAseq-tools bootstrap --target /path/to/project --force`
2) `cd /path/to/project` 후 `claude` 또는 `codex` 실행
3) `mcp__scrnaseq-tools__scan_data` 또는 `mcp__scrnaseq-tools__summarize` 사용

### 로컬 LLM로 채팅 (API 키 없음)
1) 로컬 서버 실행 (Ollama, vLLM, LM Studio 중 하나)
2) 아래 중 하나 실행:
   `scRNAseq-tools chat --provider ollama`
   `scRNAseq-tools chat --provider vllm --base-url http://localhost:8000 --model local-model`
   `scRNAseq-tools chat --provider lmstudio --base-url http://localhost:1234 --model local-model`
3) 종료하려면 `/exit`

### API 키로 채팅 (사용량 과금)
1) `export OPENAI_API_KEY=...` 또는 `export ANTHROPIC_API_KEY=...`
2) 아래 중 하나 실행:
   `scRNAseq-tools chat --provider openai --model gpt-4.1-mini`
   `scRNAseq-tools chat --provider anthropic --model claude-3-5-sonnet-20240620`
3) 종료하려면 `/exit`

### CLI 로그인으로 채팅 (Claude Code / Codex)
1) `claude` 또는 `codex` CLI가 설치 및 로그인되어 있어야 함
2) 아래 실행:
   `scRNAseq-tools chat --provider claude-cli`
   `scRNAseq-tools chat --provider codex-cli`
3) 현재 디렉토리 기준으로 `.claude` 설정을 읽음

### 로컬 REPL (LLM 없음)
1) `scRNAseq-tools repl`
2) `scan`, `summary`, `env` 사용

### 재현성 환경 (conda/Docker)
1) `conda env create -f environment.yml`
2) `conda activate scrnaseq-tools`
3) `python -m pip install -e .`

Docker:
```bash
docker build -t scrnaseq-tools:latest .
docker run --rm -it -v "$PWD:/workspace" scrnaseq-tools:latest scRNAseq-tools --help
```

## LLM 사용 (선택 사항)

기본적으로 API 호출은 없습니다. 외부 모델을 사용하려면 provider와 key를 설정하세요:

```bash
export LLM_PROVIDER=openai
export OPENAI_API_KEY=your_key_here
export LLM_MODEL=gpt-4.1-mini
```

지원 provider (기본 내장): `none`, `openai`, `anthropic`, `ollama`.

### LLM 설정 조합

Provider 결정 우선순위:
1) `LLM_PROVIDER`
2) `OPENAI_API_KEY` (자동으로 `openai`로 인식)
3) `ANTHROPIC_API_KEY` (자동으로 `anthropic`로 인식)
4) 기본값 `none`

Model 결정 우선순위:
1) `LLM_MODEL`
2) `MODEL`
3) `OPENAI_MODEL`
4) `ANTHROPIC_MODEL`

대표 조합:
- `LLM_PROVIDER=none` (또는 미설정) + API key 없음: LLM 비활성화.
- `OPENAI_API_KEY=...`만 설정: provider 자동으로 `openai`.
- `ANTHROPIC_API_KEY=...`만 설정: provider 자동으로 `anthropic`.
- `LLM_PROVIDER=ollama`: 로컬 Ollama 사용, API key 불필요.
- `LLM_PROVIDER=openai`인데 `OPENAI_API_KEY` 없음: provider는 설정되나 key는 누락.
- `LLM_PROVIDER=anthropic`인데 `ANTHROPIC_API_KEY` 없음: provider는 설정되나 key는 누락.

## 프로젝트 부트스트랩

기본 `.claude` 번들을 프로젝트에 설치:

```bash
scRNAseq-tools init --target /path/to/project
```

기본 소스 경로가 없으면 패키지에 포함된 `.claude` 템플릿을 사용합니다.

프로젝트 초기화와 MCP 서버 의존성 설치를 한 번에:

```bash
scRNAseq-tools bootstrap --target /path/to/project
```

`bootstrap`은 `pip`를 사용하므로 MCP 의존성 설치에 네트워크가 필요할 수 있습니다.

소스 디렉토리 직접 지정:

```bash
scRNAseq-tools init --source /path/to/.claude --target /path/to/project
```

`SCRNASEQ_TOOLS_CLAUDE_SOURCE` 환경변수로 기본 소스 경로를 변경할 수 있습니다.

## MCP 연동

내장된 `.claude` 템플릿에는 `scrnaseq-tools` MCP 서버가 포함되어 있으며,
`list_dir`, `scan_data`, `summarize`, `env` 도구를 제공합니다.

```bash
scRNAseq-tools bootstrap --target /path/to/project
```

이후 해당 프로젝트에서 Claude Code/Codex를 열면 `mcp__scrnaseq-tools__*` 도구를 사용할 수 있습니다.

### CLI 명령 조합

일반 명령:
```bash
scRNAseq-tools --help
scRNAseq-tools info
scRNAseq-tools doctor
scRNAseq-tools chat
scRNAseq-tools repl
scRNAseq-tools version
```

`init` 옵션과 동작:
- `--target DIR`: `DIR/.claude`에 설치 (기본값: 현재 디렉토리)
- `--source PATH`: 소스 `.claude` 디렉토리 지정
- `--force`: 기존 `DIR/.claude` 덮어쓰기
- `--source` 미지정 시 소스 결정 우선순위:
  1) `SCRNASEQ_TOOLS_CLAUDE_SOURCE`
  2) `/home/kwy7605/data_61/Vaccine_V2/.claude` (존재할 경우)
  3) 패키지에 포함된 템플릿
- `DIR/.claude`가 이미 있고 `--force`가 없으면 `init`은 에러로 종료

예시:
```bash
scRNAseq-tools init --target /path/to/project
scRNAseq-tools init --target /path/to/project --force
scRNAseq-tools init --source /path/to/.claude --target /path/to/project
SCRNASEQ_TOOLS_CLAUDE_SOURCE=/path/to/.claude scRNAseq-tools init --target /path/to/project
```

`bootstrap` 옵션과 동작:
- `--target DIR`: `init`과 동일한 대상 규칙
- `--source PATH`: `init`과 동일한 소스 규칙
- `--force`: 기존 `.claude`를 덮어쓴 뒤 MCP 의존성 설치
- `--no-mcp-install`: MCP 의존성 설치 생략
- `DIR/.claude`가 있고 `--force`가 없으면 기존 `.claude`를 사용

예시:
```bash
scRNAseq-tools bootstrap --target /path/to/project
scRNAseq-tools bootstrap --target /path/to/project --force
scRNAseq-tools bootstrap --target /path/to/project --no-mcp-install
scRNAseq-tools bootstrap --source /path/to/.claude --target /path/to/project
```

## 채팅 모드

`chat`은 먼저 provider를 선택한 뒤 대화형 세션을 시작합니다. Enter를 누르면 기본값(`claude-cli`)이 사용됩니다:

```bash
scRNAseq-tools chat
```

`--provider`로 프롬프트를 건너뛸 수 있습니다:

```bash
scRNAseq-tools chat --provider ollama
scRNAseq-tools chat --provider vllm --base-url http://localhost:8000 --model local-model
scRNAseq-tools chat --provider lmstudio --base-url http://localhost:1234 --model local-model
scRNAseq-tools chat --provider openai --model gpt-4.1-mini
scRNAseq-tools chat --provider anthropic --model claude-3-5-sonnet-20240620
```

로컬 provider:
- `ollama`: `OLLAMA_HOST` (기본 `http://localhost:11434`), `OLLAMA_MODEL` 사용
- `vllm`: OpenAI 호환; `VLLM_BASE_URL` (기본 `http://localhost:8000`), `VLLM_MODEL`
- `lmstudio`: OpenAI 호환; `LMSTUDIO_BASE_URL` (기본 `http://localhost:1234`), `LMSTUDIO_MODEL`

클라우드 provider:
- `openai`: `OPENAI_API_KEY` 필요, `OPENAI_MODEL` 또는 `LLM_MODEL` 사용
- `anthropic`: `ANTHROPIC_API_KEY` 필요, `ANTHROPIC_MODEL` 또는 `LLM_MODEL` 사용

CLI 로그인 재사용:
- `codex-cli`: `codex` 실행 (`SCRNASEQ_CODEX_CMD` 또는 `--cli-cmd`로 변경)
- `claude-cli`: `claude` 실행 (`SCRNASEQ_CLAUDE_CMD` 또는 `--cli-cmd`로 변경)
- `--cli-args`로 CLI 추가 인자를 전달 가능
- `claude`가 Codex로 감지되면 `.claude`가 있을 때 `AGENTS.md`를 자동 생성
- 자동 감지 비활성화: `SCRNASEQ_DISABLE_CLI_DETECT=1`
- AGENTS 자동 생성 비활성화: `SCRNASEQ_DISABLE_AGENTS_BOOTSTRAP=1`
- `claude-cli` 실행 시 `.claude`가 있으면 `CLAUDE.md`를 자동 생성
- CLAUDE.md 자동 생성 비활성화: `SCRNASEQ_DISABLE_CLAUDE_BOOTSTRAP=1`

채팅 명령:
- `/exit` 또는 `/quit`: 종료
- `/reset`: 히스토리 초기화
- `/help`: 명령 도움말

## 로컬 REPL

`repl`은 외부 LLM 없이 로컬 분석 콘솔을 실행합니다:

```bash
scRNAseq-tools repl
```

명령:
- `help`: 도움말
- `pwd`: 현재 디렉토리
- `cd <path>`: 디렉토리 변경
- `ls [path]`: 디렉토리 목록
- `scan [path] [-r]`: 데이터 파일 검색 (`.h5ad`, `.mtx`, `.csv` 등)
- `summary <path>`: 간단 요약 (`.h5ad`는 `anndata`, 표는 `pandas`)
- `env`: 설치된 분석 패키지 확인

## 재현성

기본적으로 모든 CLI 실행은 현재 디렉토리에 `run_manifest.json`을 생성합니다.
비활성화하거나 경로를 변경할 수 있습니다:

```bash
export SCRNASEQ_MANIFEST=0
export SCRNASEQ_MANIFEST_PATH=/path/to/run_manifest.json
```

전체 가이드는 `reproducibility.md`를 참고하세요.

## 새 명령 추가

`src/scrnaseq_tools/commands/`에 모듈을 만들고
`src/scrnaseq_tools/commands/__init__.py`에 등록하세요.

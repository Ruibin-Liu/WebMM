# Repository Agent Rules

This AGENTS.md extends the GLOBAL AGENTS.md.
All rules in this file are mandatory for this repository and may only
add stricter constraints, never weaken global rules.

Failure to follow these rules is an error.

---

## Stack Overview

Backend:
- FastAPI (async)
- Python 3.11+
- Async SQLAlchemy
- PostgreSQL
- JWT-based authentication
- Dependency and execution management via **uv**

Frontend:
- React
- TypeScript
- Material UI
- Zustand
- React Router
- Axios (centralized client with interceptors)

---

## Mandatory Workflow (In Addition to Global Rules)

For every non-trivial task:

1. Read `CODE_STATUS.md`
2. Create or overwrite `PLAN.md`
3. Implement strictly according to `PLAN.md`
4. Update `CODE_STATUS.md` after completion
5. CODE_STATUS.md and PLAN.md must follow their templates exactly. Do not add sections or free-form commentary.

Do not skip steps.
If a step cannot be completed, STOP and ask.

---

## Python Dependency Management (Strict)

This repository uses **uv** as the sole Python dependency and execution manager.

### Allowed Commands
- Add dependencies: `uv add <package>`
- Remove dependencies: `uv remove <package>`
- Sync environment: `uv sync`
- Run code or tools: `uv run <command>`

### Forbidden Commands
- `pip install`, `pipx`
- `poetry`, `pipenv`, `conda`
- `uv pip install` (pip-compatibility mode)
- Manually editing dependency files unless explicitly instructed

### Rationale
- `uv add/remove/sync` are the canonical, state-aware workflow
- `uv pip install` is a compatibility interface and MUST NOT be used in managed projects

---

## Running Backend Code & Tests

All Python execution MUST go through `uv run`.

Examples:
- Start API: `uv run fastapi dev`
- Run scripts: `uv run python scripts/task.py`
- Run tests: `uv run pytest`

Do NOT assume:
- An activated virtual environment
- System Python
- Global site-packages

---

## Backend Coding Rules (FastAPI)

- All API endpoints must be `async`
- No blocking I/O in request handlers
- Use dependency injection for:
  - Database sessions
  - Authentication
- Do not bypass security dependencies
- Follow existing error-handling patterns
- Preserve existing response models unless explicitly changing the API

---

## Frontend Rules (React + TypeScript)

- TypeScript is mandatory (no `any` unless explicitly justified)
- Prefer functional components
- State management via **Zustand only**
- All API calls must go through the shared Axios client
- Do not call backend endpoints directly
- Respect existing component and folder structure

---

## API Contract Rules

- Backend response shapes must remain stable
- Frontend changes must respect existing API contracts
- Breaking API changes require:
  - Explicit mention in `PLAN.md`
  - Documentation in `CODE_STATUS.md`

---

## Scope & Safety Constraints

- No refactors unless explicitly requested
- No dependency changes without approval
- No UI redesigns unless requested
- No speculative optimizations
- No “while we’re here” changes

---

## Completion Checklist

A task is complete only if:
- Implementation matches `PLAN.md`
- `CODE_STATUS.md` is updated accurately
- Changes are summarized clearly

# /ship — Commit, push, open PR, wait for CI, merge

Automates the full ship workflow for the current git branch:
lint → commit → push → PR → CI check → squash merge → delete branch.

Works in any git repo. Detects the correct Python interpreter automatically.

---

## Step 0 — Preflight

Refresh PATH so `gh` is available:
```powershell
$env:Path = [System.Environment]::GetEnvironmentVariable("Path","Machine") + ";" + [System.Environment]::GetEnvironmentVariable("Path","User")
```

Confirm we are NOT on `main` or `master` — abort if so and tell the user to switch to a feature branch.

Run `git status` to understand what has changed.

---

## Step 1 — Detect Python interpreter

Read `.vscode/settings.json` in the repo root and extract `python.defaultInterpreterPath`.
- If found and the file exists on disk → use it.
- If not found → fall back to `python` (system default).

```powershell
$vscode = Get-Content ".vscode/settings.json" -Raw -ErrorAction SilentlyContinue | ConvertFrom-Json
$python = $vscode.'python.defaultInterpreterPath'
if (-not $python -or -not (Test-Path $python)) { $python = "python" }
```

---

## Step 2 — Lint

Run ruff using the detected interpreter:
```powershell
& $python -m ruff check . --fix
& $python -m ruff format .
```

- If errors remain after `--fix`, **stop and report** — do not commit.
- If ruff reformatted or fixed files, re-stage them before committing.

---

## Step 3 — Commit

- If the user supplied text after `/ship`, use it as the commit message verbatim.
- Otherwise inspect `git diff --staged` and `git diff` to write a concise message focused on *why* (not just what).
- Stage modified tracked files by name (never `git add .` or `git add -A`).
- Commit with PowerShell here-string:
  ```powershell
  git commit -m @'
  <message>

  Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>
  '@
  ```
- If the pre-commit hook reformats files and the commit fails, re-stage the reformatted files and commit again automatically — do this at most once.

---

## Step 4 — Push

```powershell
git push origin <branch>
```

Print the remote URL on success.

---

## Step 5 — PR

Check if a PR already exists for this branch:
```powershell
gh pr list --head <branch> --json number,url,title
```

- If a PR exists → print its URL and use its number for CI polling.
- If no PR exists → create one:
  ```powershell
  gh pr create --title "<title>" --body "<body>" --base main
  ```
  Use the branch name and commit log to write the title and body. Print the new PR URL.

---

## Step 6 — Wait for CI

Poll every 30 seconds:
```powershell
gh pr checks <number> --json name,state,conclusion
```

Print a one-line update each poll: `[HH:MM] Waiting — N/M checks complete`

Stop when ALL checks are `COMPLETED`. Then evaluate:

- **`codecov/patch`** → non-blocking, treat as warning only.
- Any other check with `conclusion = FAILURE` → print the failed check names and their `detailsUrl`, **stop, do not merge**, ask the user whether to fix and retry.
- All required checks `SUCCESS` → proceed to merge.

---

## Step 7 — Merge

```powershell
gh pr merge <number> --squash --delete-branch
```

Confirm success, print the resulting commit SHA on `main`, then:
```powershell
git fetch origin main
```

---

## Notes

- Always use the **PowerShell tool** (never Bash) and the `&` call operator for executables.
- If `gh` auth fails, instruct the user to run `gh auth login` in a new terminal and retry.
- If there are no staged or unstaged changes, skip Steps 2–4 and go straight to Step 5 (the branch may already be pushed with an open PR).

<!-- markdownlint-disable MD041 -->
## TrackTx Solo Git Workflow (serhataktay/TrackTx)

This guide describes a clear, repeatable daily workflow for developing TrackTx solo, across one or multiple computers. It includes commands, what they do, and how to handle common situations safely.

Audience: You (solo maintainer). Remote: `github.com/serhataktay/TrackTx`.

---

### 1) Prerequisites

- Install Git on each computer you use.
- Set up SSH keys with GitHub for password‑less pushes (recommended). Otherwise, use HTTPS.
- Optional: Install Git LFS if you expect large binary files (e.g., datasets, images).

Verify Git:

```bash
git --version
```

If you’ll use SSH, ensure you can reach GitHub:

```bash
ssh -T git@github.com
```

If this greets you by username, your SSH is configured.

---

### 2) One‑time Git configuration (per computer)

These settings make daily work smoother and keep history clean.

```bash
git config --global user.name "Serhat Aktay"
git config --global user.email "YOUR_EMAIL@example.com"   # Use your GitHub email

# Keep pulls linear (avoids unnecessary merge commits)
git config --global pull.rebase true

# When pushing a new local branch, automatically set upstream
git config --global push.autoSetupRemote true

# Optional: default branch name when initializing new repos
git config --global init.defaultBranch main

# Helpful aliases
git config --global alias.st "status -sb"
git config --global alias.co switch
git config --global alias.br "branch -vv"
git config --global alias.lg "log --oneline --graph --decorate --all"
```

What these do:

- `pull.rebase true`: Your local commits are replayed on top of the latest remote commits instead of creating a merge commit. Cleaner history, fewer conflicts later.
- `push.autoSetupRemote true`: First push of a new branch automatically links it to the branch on GitHub.
- Aliases: Shortcuts for common commands.

---

### 3) Cloning the repository

Do this once per computer.

SSH (recommended):

```bash
git clone git@github.com:serhataktay/TrackTx.git
cd TrackTx
```

HTTPS (if SSH not set up):

```bash
git clone https://github.com/serhataktay/TrackTx.git
cd TrackTx
```

---

### 4) Branching model

Even solo, use short‑lived branches for tasks. This gives you safety (easy rollback), organization, and clean main history.

- `main`: Always deployable/stable.
- Feature/fix branches: `feat/<task>`, `fix/<issue>`, `chore/<maintenance>`.

Examples:

- `feat/add-qc-docs`
- `fix/normalize-tracks-argparse-bug`
- `chore/update-deps`

---

### 5) Daily routine

Perform these steps at the start, during, and end of your day.

#### Start of day: sync and branch

```bash
cd /Volumes/samsung/TrackTx   # or the path where you cloned
git switch main
git pull --rebase origin main       # get latest and rebase local changes

# Create or switch to your task branch
git switch -c feat/<short-task-name>   # if starting a new task branch
# or
git switch feat/<short-task-name>      # if the branch already exists locally
```

What this does:

- Ensures `main` is up to date before you branch.
- Keeps your history linear and conflict‑free as long as possible.

#### During the day: work in small commits, push to back up

```bash
# Stage changes interactively (lets you commit logically related chunks)
git add -p

# Write a concise message describing what changed
git commit -m "feat: concise description of change"

# Publish your branch or update remote
git push
```

Notes:

- Commit saves locally. Push backs up to GitHub (and makes the work accessible from other computers).
- Prefer multiple small commits over one giant commit; they’re easier to review and revert.

#### End of day: ensure everything is on GitHub

```bash
git status   # verify nothing important is left unstaged/uncommitted
git push     # make sure all commits are on GitHub
```

---

### 6) Picking up work on a different computer

First time on that computer:

```bash
git clone git@github.com:serhataktay/TrackTx.git
cd TrackTx
git fetch                      # updates the list of remote branches
git switch feat/<short-task-name>
git pull --rebase              # pull the latest for that branch
```

If the repo already exists on that computer:

```bash
cd TrackTx
git fetch
git switch feat/<short-task-name>
git pull --rebase
```

Why this works:

- Because you pushed from computer A, computer B can fetch and switch to the same branch name and continue seamlessly.

---

### 7) Finishing a task (merge to main)

When a task is complete, tested, and you want it on `main`:

```bash
# Ensure main is up to date
git switch main
git pull --rebase origin main

# Rebase your branch on top of latest main (to ensure clean, fast-forward merge)
git switch feat/<short-task-name>
git rebase main

# Merge into main using fast-forward only (prevents accidental merge commits)
git switch main
git merge --ff-only feat/<short-task-name>
git push

# Optionally delete the branch locally and remotely
git branch -d feat/<short-task-name>
git push origin --delete feat/<short-task-name>
```

What this does:

- `rebase main`: Replay your branch commits on top of current `main` to avoid merge commits.
- `merge --ff-only`: Fails if a fast‑forward is not possible, protecting history. If it fails, re‑run the rebase step.

---

### 8) Keeping your branch up to date with main (while in progress)

If `main` changed while you are still working:

```bash
git switch main
git pull --rebase origin main
git switch feat/<short-task-name>
git rebase main

# If conflicts appear:
#   1) Manually edit conflicted files to resolve
#   2) git add <resolved-file>
#   3) git rebase --continue
```

Why rebase instead of merge here:

- Rebase keeps a straight line of history, which makes future merges cleaner and logs easier to read.

---

### 9) Resolving conflicts (quick guide)

When Git reports conflicts during `pull --rebase` or `rebase`:

1. Open each conflicted file. Git marks conflicts with `<<<<<<<`, `=======`, `>>>>>>>`.
2. Decide which changes to keep (yours, theirs, or a manual combination).
3. Remove conflict markers, ensure the code builds/tests.
4. Stage resolved files: `git add <file>`.
5. Continue: `git rebase --continue`.
6. If it gets messy: `git rebase --abort` to go back to the pre‑rebase state and ask for help or try again carefully.

---

### 10) Recovering from mistakes safely

- Show recent history: `git lg` (alias) or `git log --oneline --graph --decorate --all`.
- Undo last commit but keep changes staged: `git reset --soft HEAD~1`.
- Undo last commit and unstage changes: `git reset --mixed HEAD~1`.
- Discard local uncommitted changes in a file: `git checkout -- <file>` (dangerous; cannot undo).
- Abort a problematic rebase: `git rebase --abort`.

Tip: When in doubt, create a safety branch before risky operations:

```bash
git switch -c safety/backup-before-rebase
```

---

### 11) Tags and releases (optional but useful)

When you reach a stable milestone, create an annotated tag:

```bash
git switch main
git pull --rebase origin main
git tag -a vX.Y.Z -m "Description of this release"
git push origin vX.Y.Z
```

Tags make it easy to reproduce runs or share specific stable states.

---

### 12) Large files (Git LFS)

If you need to track large binary files (e.g., `.bam`, `.bigWig`, `.fastq` beyond GitHub’s standard size limits), configure Git LFS once:

```bash
git lfs install
git lfs track "*.bigWig"
git add .gitattributes
git commit -m "chore: track .bigWig with Git LFS"
git push
```

Then add/push large files as usual. LFS handles storage efficiently.

---

### 13) SSH vs HTTPS

- SSH: `git@github.com:serhataktay/TrackTx.git` — no password after initial key setup.
- HTTPS: `https://github.com/serhataktay/TrackTx.git` — may prompt for credentials/tokens.

Prefer SSH for convenience and reliability.

---

### 14) Common commands cheat‑sheet

```bash
# Status and branches
git st                         # alias: status -sb
git br                         # alias: branch -vv
git lg                         # alias: log graph

# Sync main
git switch main
git pull --rebase origin main

# Start/continue a task
git switch -c feat/<task>      # new branch
git switch feat/<task>         # existing branch
git add -p
git commit -m "feat: ..."
git push

# Update branch with latest main
git switch main && git pull --rebase
git switch feat/<task>
git rebase main

# Finish task
git switch main
git merge --ff-only feat/<task>
git push
git branch -d feat/<task>
git push origin --delete feat/<task>

# Conflict recovery
git rebase --continue
git rebase --abort
```

---

### 15) Daily checklist

- Start of day:
  - Pull latest `main`.
  - Create/switch to your task branch.
- During the day:
  - Commit in small, logical chunks.
  - Push regularly to back up work.
- End of day:
  - Ensure `git status` is clean or intentionally dirty.
  - `git push` so other computers (and future you) can pick up.

---

### 16) FAQs

**Q: I committed but don’t see changes on GitHub.**
A: You need to `git push` to publish local commits.

**Q: I want to move to another computer. What’s required?**
A: Make sure your branch is pushed (`git push`). On the other computer, `git fetch`, `git switch <branch>`, `git pull --rebase`.

**Q: Rebase showed conflicts; I’m stuck.**
A: Resolve files, `git add <file>`, then `git rebase --continue`. If overwhelmed, `git rebase --abort` and try again in smaller steps.

**Q: Should I ever commit on `main` directly?**
A: Prefer not to. Use branches and fast‑forward merges to keep `main` clean and stable.

---

With this workflow, your work stays organized, easy to synchronize across machines, and resilient to mistakes.

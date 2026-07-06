# Repo conventions — hydroSense

## Git / merge policy

- **Do not squash-merge PRs.** Merge with a merge commit
  (`gh pr merge --merge`) or rebase (`gh pr merge --rebase`). Squashing
  collapses a branch’s commits into a single new commit on `main`, so
  `git branch --merged main` can no longer tell the branch is merged —
  every squash-merged branch then shows as “unmerged / ahead” and piles
  up as stale clutter that is tedious and risky to audit later. Squash
  merging is **disabled in the repo settings**; keep it that way.

- **Delete each branch as soon as its PR is merged** (remote *and*
  local). The repo has *Automatically delete head branches* enabled, so
  merging a PR removes the remote branch; delete your local copy too
  (`git branch -d <branch>`). Don’t leave merged branches lying around.

- Merge commits and rebase merges are both allowed; pick whichever gives
  the cleaner history for the change.

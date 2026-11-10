# Create a separate GitHub repository

This archive is a **new standalone repository source tree**; it has **not** been pushed to GitHub. Do **not** copy the parent project's `.git/` directory.

1. Create a new empty GitHub repository, for example `late-time-tails-5plus1-mmief`, and leave README/license creation unchecked (this project already includes them).
2. Extract the ZIP and open the `late_time_tails_5plus1_mmief` folder in VS Code or PowerShell.
3. Run:

```powershell
git init
git add .
git commit -m "Initial six-dimensional MMIEF research prototype"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/late-time-tails-5plus1-mmief.git
git push -u origin main
```

Replace `YOUR_USERNAME` with your GitHub username. Run the local Julia tests before publishing scientific results. The repository deliberately ignores generated output under `results/`.

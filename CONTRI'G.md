# How to Upload Your Code to the BENG207 GitHub Repository

## A Step-by-Step Guide for BENG207 students

This guide will walk you through uploading your C code files to the shared class repository on GitHub.

---

## Prerequisites

Before you begin, make sure you have:

1. A **GitHub account** (sign up free at [github.com](https://github.com))
2. **Git** installed on your computer
3. Your C code files ready to upload

### Check if Git is Installed

Open your terminal (Windows is okay) and type:

```bash
git --version
```

If you see a version number, you're good to go. If not, install Git:

- **Fedora/RHEL:** `sudo dnf install git`
- **Ubuntu/Debian:** `sudo apt install git`
- **macOS:** `xcode-select --install`
- **Windows:** Download from [git-scm.com](https://git-scm.com)

---

## Step 1: Configure Git (First Time Only)

Set your name and email so your contributions are properly attributed:

```bash
git config --global user.name "Your Full Name (Anna ...`)"
git config --global user.email "your_email@ucsd.edu"
```

---

## Step 2: Create a Personal Access Token on GitHub

GitHub no longer accepts passwords for Git operations. You need a Personal Access Token (PAT):

1. Go to **github.com** and log in
2. Click your **profile picture** (top right) → **Settings**
3. Scroll down the left sidebar → **Developer settings**
4. Click **Personal access tokens** → **Tokens (classic)**
5. Click **Generate new token** → **Generate new token (classic)**
6. Fill in the details:
   - **Note:** Give it a name like `my-laptop-git`
   - **Expiration:** Choose an expiration (e.g., 90 days)
   - **Scopes:** Check the **`repo`** box (full control of repositories)
7. Click **Generate token**
8. **Copy the token immediately** — you will not be able to see it again!

Save this token somewhere safe. You will use it as your "password" when pushing code.

---

## Step 3: Fork the Repository (Recommended Method)

Forking creates your own copy of the repository where you can freely make changes.

1. Go to: [https://github.com/akmatsuoka/BENG207_ORGANOIDS_MODEL](https://github.com/akmatsuoka/BENG207_ORGANOIDS_MODEL)
2. Click the **Fork** button (top right of the page)
3. This creates a copy under your own GitHub account

---

## Step 4: Clone Your Fork to Your Computer

```bash
cd ~/Projects
git clone https://github.com/YOUR_USERNAME/BENG207_ORGANOIDS_MODEL.git
cd BENG207_ORGANOIDS_MODEL
```

Replace `YOUR_USERNAME` with your actual GitHub username.

When prompted:

- **Username:** Your GitHub username
- **Password:** Paste your **Personal Access Token** (not your GitHub password)

---

## Step 5: Add Your Code Files

Copy your C files into the cloned repository folder. You can do this with the file manager or from the terminal:

```bash
cp ~/path/to/your/code_files/*.c ~/Projects/BENG207_ORGANOIDS_MODEL/
cp ~/path/to/your/code_files/*.h ~/Projects/BENG207_ORGANOIDS_MODEL/
```

Adjust the paths to match where your files are located.

---

## Step 6: Stage, Commit, and Push Your Changes

Run these commands one at a time:

```bash
cd ~/Projects/BENG207_ORGANOIDS_MODEL

git add .

git commit -m "Add my organoid model code - [Your Name]"

git push origin main
```

When prompted:

- **Username:** Your GitHub username
- **Password:** Paste your **Personal Access Token**

---

## Step 7: Create a Pull Request

A pull request lets the instructor review your code before merging it into the main repository.

1. Go to your fork on GitHub: `https://github.com/YOUR_USERNAME/BENG207_ORGANOIDS_MODEL`
2. Click **Contribute** → **Open pull request**
3. Add a title describing your changes (e.g., "Add organoid simulation code - Jane Doe")
4. Add any notes or comments in the description
5. Click **Create pull request**

Dr. Matsuoka or Kevin will review your code and merge your code.

---

## Quick Reference: Common Commands

| What you want to do        | Command                        |
| -------------------------- | ------------------------------ |
| Check status of your files | `git status`                   |
| See what changed           | `git diff`                     |
| Stage all changes          | `git add .`                    |
| Stage a specific file      | `git add filename.c`           |
| Commit with a message      | `git commit -m "your message"` |
| Push to GitHub             | `git push origin main`         |
| Pull latest changes        | `git pull origin main`         |

---

## Troubleshooting

**"Repository not found" error:**
Make sure you have forked the repository and are pushing to your fork, not the instructor's repository.

**"Authentication failed" error:**
You need to use a Personal Access Token, not your GitHub password. See Step 2 above.

**"rejected - fetch first" error:**
Run `git pull origin main` first, then try pushing again.

**"nothing to commit" message:**
Your files are either already committed or not in the repository folder. Check with `git status`.

---

## Save Your Credentials (Optional)

To avoid entering your token every time:

```bash
git config --global credential.helper store
```

Then push once — your credentials will be saved for future use.

---

## Need Help?

Contact your instructor or visit the [GitHub Docs](https://docs.github.com) for more detailed guides or comtact Kevin Nella or Dr. Matsuoka directly. 

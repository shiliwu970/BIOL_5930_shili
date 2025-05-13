# Install the `usethis` and `gitcreds` packages
install.packages(c("usethis", "gitcreds"))
library(usethis); library(gitcreds)

# Add your GitHub username and email
usethis::use_git_config(user.name = "shiliwu790",
                        user.email = "shili.wu@slu.edu")

# Create a token (Note this will open a GitHub browser tab)
## See steps 6-10 in GitHub's tutorial (link below)
usethis::create_github_token()

# Now, give your token to RStudio
gitcreds::gitcreds_set()
## After you run this line you'll follow some prompts in the "Console" tab of RStudio
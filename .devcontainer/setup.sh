#!/usr/bin/env bash
# Devcontainer post-create setup. Runs once after Codespaces builds the
# container. Installs the two coding-agent CLIs participants use during the
# workshop and pre-configures Pi to default to OpenRouter.
set -euo pipefail

# Claude Code CLI. Routed to OpenRouter via ANTHROPIC_AUTH_TOKEN (set in
# devcontainer.json's remoteEnv from the OPENROUTER_API_KEY Codespaces secret).
npm install -g @anthropic-ai/claude-code

# Pi coding agent. Reads OPENROUTER_API_KEY directly from the environment.
# --ignore-scripts disables npm lifecycle scripts on install, per the
# maintainer's published install instructions at pi.dev.
npm install -g --ignore-scripts @earendil-works/pi-coding-agent

# Default Pi to OpenRouter + DeepSeek V4 Flash so participants don't have to
# /login or pick a provider on first launch, and so the cheap default model
# is in place from day one. Pi writes back to this file when the user
# switches providers/models via Ctrl+L or /model, so anything participants
# change persists across the session.
mkdir -p "$HOME/.pi/agent"
cat > "$HOME/.pi/agent/settings.json" <<'JSON'
{
  "defaultProvider": "openrouter",
  "defaultModel": "deepseek/deepseek-v4-flash"
}
JSON

# Pi plugins for the workshop. Installed via `pi install`, which writes to
# $HOME/.pi/agent/packages/. All five work alongside Pi's built-in read /
# write / edit / bash tools, none of which require plugins.
#
#  - pi-web-access: web search, URL fetching, PDF extraction. For looking up
#    R / tidyverse / Quarto docs from inside Pi.
#  - rpiv-ask-user-question: structured questions instead of free-form
#    guessing when context is ambiguous. Beginner-friendly.
#  - rpiv-todo: persistent todo overlay that survives /reload and
#    conversation compaction. Fits the workshop's three-part exercise shape.
#  - rpiv-advisor: lets Pi request a second opinion from a stronger
#    reviewer model before acting. Demonstrates the cost/strength tradeoff
#    explicitly (cheap default + smart-model escalation), which is one of
#    the workshop's threads.
#  - rpiv-i18n: adds the /languages command so participants whose first
#    language isn't English can switch the rpiv-* plugin UIs to French
#    (or other supported locales).
pi install npm:pi-web-access
pi install npm:@juicesharp/rpiv-ask-user-question
pi install npm:@juicesharp/rpiv-todo
pi install npm:@juicesharp/rpiv-advisor
pi install npm:@juicesharp/rpiv-i18n

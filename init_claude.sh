#!/bin/bash
# Claude Code 环境初始化（bash 版，供 Hermes 会话调用；与 init_claude.ps1 同源）
export ANTHROPIC_BASE_URL="https://open.bigmodel.cn/api/anthropic"
export ANTHROPIC_AUTH_TOKEN="4644e4d94ccc4b09a0de5dc3036b6e4f.s5f7JH5sCegjWGu1"
export ANTHROPIC_MODEL="glm-5.3"
export ANTHROPIC_DEFAULT_OPUS_MODEL="glm-5.3"
export ANTHROPIC_DEFAULT_SONNET_MODEL="glm-5.3"
export ANTHROPIC_DEFAULT_HAIKU_MODEL="glm-5-turbo"
export CLAUDE_CODE_SUBAGENT_MODEL="glm-5-turbo"
export CLAUDE_CODE_DISABLE_UNKNOWN_MODEL_WINDOW_ENFORCEMENT=1

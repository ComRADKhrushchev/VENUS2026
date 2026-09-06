<#
.SYNOPSIS
Claude Code 启动脚本 (PowerShell 版本)
.DESCRIPTION
设置环境变量并启动 claude --bare，适配 Deepseek 兼容端点。
请修改下面的 $apiKey 和 $baseUrl 变量。
#>

# --- 配置区域，请修改为你自己的值 ---
$apiKey  = "sk-c50d68aaee464eaa9fecbf3d9435467c"
$baseUrl = "https://api.deepseek.com/anthropic"
# 可选：指定模型
$env:ANTHROPIC_BASE_URL="https://open.bigmodel.cn/api/anthropic"
$env:ANTHROPIC_AUTH_TOKEN="4644e4d94ccc4b09a0de5dc3036b6e4f.s5f7JH5sCegjWGu1"
$env:ANTHROPIC_MODEL="glm-5.3"
$env:ANTHROPIC_DEFAULT_OPUS_MODEL="glm-5.3"
$env:ANTHROPIC_DEFAULT_SONNET_MODEL="glm-5.3"
$env:ANTHROPIC_DEFAULT_HAIKU_MODEL="glm-5-turbo"
$env:CLAUDE_CODE_SUBAGENT_MODEL="glm-5-turbo"
$env:CLAUDE_CODE_EFFORT_LEVEL="max"

# --- 设置环境变量 ---
$env:ANTHROPIC_API_KEY  = $apiKey
$env:ANTHROPIC_BASE_URL = $baseUrl
if ($model) { $env:ANTHROPIC_DEFAULT_MODEL = $model }

# --- 启动 Claude（bare 模式 + 清空系统提示词）---
claude 
#!/usr/bin/env bash
# Script to pull Nextflow pipeline information from ebi-metagenomics organization
# Requires: gh (GitHub CLI) installed and configured

set -euo pipefail

ORG="EBI-Metagenomics"

echo "Fetching all repositories from $ORG organization..."
echo ""

# Get all repositories (excluding nf-modules)
# Using a high limit since gh repo list doesn't support proper pagination with --page
all_repos=$(gh repo list "$ORG" --limit 500 --json name --jq '.[].name' 2>/dev/null || echo "")

# Filter out nf-modules and empty lines
all_repos=$(echo "$all_repos" | grep -v '^nf-modules$' | grep -v '^$' || echo "")

if [ -z "$all_repos" ]; then
  echo "No repositories found"
  exit 0
fi

echo "Checking $(echo "$all_repos" | wc -l | tr -d ' ') repositories for Nextflow pipelines..."
echo ""

# Temporary file to store pipeline data as JSON
temp_file=$(mktemp)

for repo in $all_repos; do
  # Get repository info
  repo_info=$(gh api "repos/$ORG/$repo" 2>/dev/null || echo "")

  if [ -z "$repo_info" ]; then
    continue
  fi

  # Check if repository is private
  is_private=$(echo "$repo_info" | jq -r '.private')
  if [ "$is_private" = "true" ]; then
    continue
  fi

  # Get repository root contents
  repo_contents=$(gh api "repos/$ORG/$repo/contents" --jq '.[].name' 2>/dev/null || echo "")

  # Check if it's a Nextflow pipeline (has main.nf or nextflow.config)
  has_nextflow=$(echo "$repo_contents" | grep -E '^(main\.nf|nextflow\.config)$' || echo "")

  if [ -z "$has_nextflow" ]; then
    continue
  fi

  # Extract info
  description=$(echo "$repo_info" | jq -r '.description // "No description"')
  stars=$(echo "$repo_info" | jq -r '.stargazers_count')
  default_branch=$(echo "$repo_info" | jq -r '.default_branch')

  # Check for modules.json
  has_modules_json=$(echo "$repo_contents" | grep -E '^modules\.json$' || echo "")

  if [ -n "$has_modules_json" ]; then
    has_modules="true"
  else
    has_modules="false"
  fi

  # Output as JSON line
  echo "{\"name\":\"$repo\",\"repo\":\"$ORG/$repo\",\"description\":\"$description\",\"stars\":$stars,\"branch\":\"$default_branch\",\"has_modules_json\":$has_modules}" >> "$temp_file"

  echo "=== $repo ==="
  echo "  Description: $description"
  echo "  Stars: $stars"
  echo "  Has modules.json: $has_modules"
  echo ""
done

echo ""
echo "Generating pipelines.yml..."
echo ""

# Sort pipelines: first by has_modules_json (true first), then by stars (descending)
sorted_pipelines=$(cat "$temp_file" | jq -s 'sort_by(.has_modules_json == false, -.stars)')

# Get script directory to write to correct location
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
OUTPUT_FILE="$SCRIPT_DIR/pipelines.yml"

# Generate YAML output using jq
{
  echo "pipelines:"
  echo "$sorted_pipelines" | jq -r '.[] | "  - name: \(.name)\n    repo: \(.repo)\n    branch: \(.branch)\n    description: \"\(.description)\"\n    has_modules_json: \(.has_modules_json)"'
} > "$OUTPUT_FILE"

rm "$temp_file"

echo ""
echo "Found $(echo "$sorted_pipelines" | jq 'length') Nextflow pipeline(s) total"
echo "Updated: $OUTPUT_FILE"
echo "Done!"

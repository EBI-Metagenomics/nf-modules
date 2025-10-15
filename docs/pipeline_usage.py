#!/usr/bin/env python3
import json
from pathlib import Path
from typing import Any

import httpx
import yaml


class PipelineUsageFetcher:
    """Fetch and parse pipeline usage information from modules.json files."""

    def __init__(self, config_path: str):
        """Initialize the fetcher.

        :param config_path: Path to pipelines.yml configuration file
        :type config_path: str
        """
        self.config_path = Path(config_path)
        self.pipelines = self._load_config()

    def _load_config(self) -> list[dict[str, str]]:
        """Load pipeline configuration from YAML file.

        :returns: List of pipeline configurations
        :rtype: list[dict[str, str]]
        """
        if not self.config_path.exists():
            print(f"Warning: Pipeline config not found at {self.config_path}")
            return []

        with open(self.config_path) as f:
            config = yaml.safe_load(f)
            return config.get("pipelines", [])

    def _fetch_modules_json(self, repo: str, branch: str) -> dict[str, Any] | None:
        """Fetch modules.json from a GitHub repository.

        :param repo: Repository in format 'owner/repo'
        :type repo: str
        :param branch: Branch name
        :type branch: str
        :returns: Parsed modules.json or None if fetch fails
        :rtype: dict[str, Any] | None
        """
        url = (
            f"https://raw.githubusercontent.com/{repo}/refs/heads/{branch}/modules.json"
        )

        try:
            with httpx.Client(timeout=10.0) as client:
                response = client.get(url)
                response.raise_for_status()
                return response.json()
        except (httpx.HTTPError, json.JSONDecodeError) as e:
            print(f"Warning: Failed to fetch modules.json from {repo} ({branch}): {e}")
            return None

    def build_module_usage_index(self) -> dict[str, list[dict[str, str]]]:
        """Build reverse index of which pipelines use which modules.

        :returns: Dictionary mapping module names to list of pipeline info
        :rtype: dict[str, list[dict[str, str]]]
        """
        usage_index = {}

        for pipeline in self.pipelines:
            name = pipeline.get("name")
            repo = pipeline.get("repo")
            branch = pipeline.get("branch", "main")
            description = pipeline.get("description", "")

            print(f"Fetching module usage from {name}...")

            modules_json = self._fetch_modules_json(repo, branch)
            if not modules_json:
                continue

            # Find the ebi-metagenomics repo key (can be different URL formats)
            # Supports: https://github.com/EBI-Metagenomics/nf-modules
            #           https://github.com/EBI-Metagenomics/nf-modules.git
            #           git@github.com:EBI-Metagenomics/nf-modules.git
            ebi_repo_key = None
            for repo_url in modules_json.get("repos", {}).keys():
                repo_url_lower = repo_url.lower()
                if (
                    "ebi-metagenomics" in repo_url_lower
                    and "nf-modules" in repo_url_lower
                ):
                    ebi_repo_key = repo_url
                    break

            if not ebi_repo_key:
                print(f"  No ebi-metagenomics/nf-modules found in {name}")
                continue

            repo_data = modules_json.get("repos", {}).get(ebi_repo_key, {})

            # Extract ebi-metagenomics modules
            ebi_modules = repo_data.get("modules", {}).get("ebi-metagenomics", {})

            for module_name in ebi_modules.keys():
                if module_name not in usage_index:
                    usage_index[module_name] = []

                usage_index[module_name].append(
                    {
                        "name": name,
                        "repo": repo,
                        "description": description,
                    }
                )

            # Extract ebi-metagenomics subworkflows
            ebi_subworkflows = repo_data.get("subworkflows", {}).get(
                "ebi-metagenomics", {}
            )

            for subworkflow_name in ebi_subworkflows.keys():
                if subworkflow_name not in usage_index:
                    usage_index[subworkflow_name] = []

                usage_index[subworkflow_name].append(
                    {
                        "name": name,
                        "repo": repo,
                        "description": description,
                    }
                )

        return usage_index

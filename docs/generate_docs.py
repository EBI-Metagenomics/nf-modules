#!/usr/bin/env python3
import shutil
from pathlib import Path
from typing import Any

import click
import yaml
from jinja2 import Environment, FileSystemLoader, Template
from pipeline_usage import PipelineUsageFetcher


class ModuleParser:
    """Parse nf-core module meta.yml files and generate documentation."""

    modules_repo_path: Path
    output_dir: Path
    pipeline_usage_index: dict[str, list[dict[str, str]]]
    jinja_env: Environment

    def __init__(
        self,
        modules_repo_path: str,
        output_dir: str,
        pipeline_usage_index: dict[str, list[dict[str, str]]] | None = None,
    ) -> None:
        """Initialize the parser.

        :param modules_repo_path: Path to the nf-core modules repository
        :type modules_repo_path: str
        :param output_dir: Output directory for generated documentation
        :type output_dir: str
        :param pipeline_usage_index: Dictionary mapping module names to pipelines
        :type pipeline_usage_index: dict[str, list[dict[str, str]]] | None
        """
        self.modules_repo_path = Path(modules_repo_path)
        self.output_dir = Path(output_dir)
        self.pipeline_usage_index = pipeline_usage_index or {}
        self.jinja_env = Environment(
            loader=FileSystemLoader("templates"), trim_blocks=True, lstrip_blocks=True
        )

    def find_meta_files(self, directory: str) -> list[Path]:
        """Find all meta.yml files in a directory.

        :param directory: Directory to search (modules or subworkflows)
        :type directory: str
        :returns: List of meta.yml file paths
        :rtype: list[Path]
        """
        search_path: Path = self.modules_repo_path / directory
        return list(search_path.glob("**/meta.yml"))

    def parse_meta_yml(self, meta_file: Path) -> dict[str, Any]:
        """Parse a meta.yml file.

        :param meta_file: Path to meta.yml file
        :type meta_file: Path
        :returns: Parsed YAML content as dictionary
        :rtype: dict[str, Any]
        """
        with open(meta_file) as f:
            return yaml.safe_load(f)

    def get_module_path_info(self, meta_file: Path, base_type: str) -> dict[str, str]:
        """Extract path information from meta.yml file location.

        :param meta_file: Path to meta.yml file
        :type meta_file: Path
        :param base_type: 'modules' or 'subworkflows'
        :type base_type: str
        :returns: Dictionary with path components
        :rtype: dict[str, str]
        """
        parts: tuple[str, ...] = meta_file.relative_to(
            self.modules_repo_path / base_type
        ).parts[:-1]
        return {
            "full_path": "/".join(parts),
            # Name used by "nf-core modules/subworkflows install"
            "component": "/".join(parts[1:]),
            "category": parts[0] if parts else "",
            "name": parts[-1] if parts else "",
            "subcategory": "/".join(parts[1:-1]) if len(parts) > 2 else "",
        }

    def _detect_component_types(
        self, components: list[str | dict[str, Any]]
    ) -> list[dict[str, str]]:
        """Detect whether components are modules or subworkflows.

        :param components: List of component names or dicts with component info
        :type components: list[str | dict[str, Any]]
        :returns: List of component info with type detection
        :rtype: list[dict[str, str]]
        """
        component_info: list[dict[str, str]] = []

        for component in components:
            # Handle both string format and dict format (for third-party modules)
            component_name: str | None
            git_remote: str | None

            if isinstance(component, dict):
                # Component is a dict like {"module_name": {"git_remote": "..."}}
                component_name = next(iter(component), None)
                if not component_name:
                    continue
                options = component[component_name]
                git_remote = (
                    options.get("git_remote") if isinstance(options, dict) else None
                )
            else:
                # Component is a simple string
                component_name = component
                git_remote = None

            # Check if component exists as a subworkflow
            subworkflow_path: Path = (
                self.modules_repo_path
                / "subworkflows"
                / "ebi-metagenomics"
                / component_name
                / "meta.yml"
            )

            component_type: str
            if subworkflow_path.exists():
                component_type = "subworkflow"
            else:
                component_type = "module"

            comp_info: dict[str, str] = {"name": component_name, "type": component_type}
            if git_remote:
                comp_info["git_remote"] = git_remote
                if "nf-core/modules" in git_remote:
                    comp_info["url"] = (
                        f"https://nf-co.re/modules/{component_name.replace('/', '_')}"
                    )

            component_info.append(comp_info)

        return component_info

    def _flatten_elements(self, value: Any) -> list[tuple[str, dict[str, Any]]]:
        """Flatten nested lists of single-key dicts into (name, info) pairs.

        nf-core meta.yml files describe a channel as a list (or a list of tuples)
        of dicts, for example ``[[{meta: {...}}, {"*.tsv": {...}}]]``.

        :param value: Input or output section, or one channel of it
        :type value: Any
        :returns: List of (element name, element info) pairs
        :rtype: list[tuple[str, dict[str, Any]]]
        """
        if isinstance(value, list):
            return [pair for item in value for pair in self._flatten_elements(item)]
        if isinstance(value, dict):
            return [(k, v) for k, v in value.items() if isinstance(v, dict)]
        return []

    def _summarise_channel(self, value: Any) -> dict[str, Any]:
        """Merge the elements of one output channel into one table row.

        The meta map is skipped when the channel has other elements.

        :param value: Channel definition from meta.yml
        :type value: Any
        :returns: Dictionary with type, description and pattern
        :rtype: dict[str, Any]
        """
        # Old flat format: {"channel": {"type": ..., "description": ...}}
        if isinstance(value, dict) and ("type" in value or "description" in value):
            return value

        elements = self._flatten_elements(value)
        infos = [info for name, info in elements if not name.startswith("meta")]
        infos = infos or [info for _, info in elements]
        if not infos:
            return {}

        return {
            "type": infos[-1].get("type", "file"),
            "description": "; ".join(
                str(info["description"]).strip()
                for info in infos
                if info.get("description")
            ),
            "pattern": infos[-1].get("pattern"),
        }

    def _process_meta_data(self, meta_data: dict[str, Any]) -> dict[str, Any]:
        """Process meta data to handle list structures in input/output.

        :param meta_data: Raw meta data from YAML
        :type meta_data: dict[str, Any]
        :returns: Processed meta data with flattened input/output structures
        :rtype: dict[str, Any]
        """
        processed: dict[str, Any] = meta_data.copy()

        # Process tools section if it's a list
        if "tools" in processed and isinstance(processed["tools"], list):
            processed_tools: dict[str, Any] = {}
            for item in processed["tools"]:
                if isinstance(item, dict):
                    processed_tools.update(item)
            for tool_info in processed_tools.values():
                # Some meta.yml files have the licence as a string, not a list
                if isinstance(tool_info, dict) and isinstance(
                    tool_info.get("licence"), str
                ):
                    tool_info["licence"] = [tool_info["licence"]]
            processed["tools"] = processed_tools

        # Process input section: one row per channel element
        if "input" in processed and isinstance(processed["input"], list):
            processed["input"] = dict(self._flatten_elements(processed["input"]))

        # Process output section: one row per channel, in list or dict format
        if "output" in processed:
            channels: list[tuple[str, Any]] = []
            if isinstance(processed["output"], list):
                for item in processed["output"]:
                    if isinstance(item, dict):
                        channels.extend(item.items())
            elif isinstance(processed["output"], dict):
                channels = list(processed["output"].items())
            processed["output"] = {
                name: self._summarise_channel(value) for name, value in channels
            }

        # Process components to add type information
        if "components" in processed and isinstance(processed["components"], list):
            processed["components"] = self._detect_component_types(
                processed["components"]
            )

        return processed

    def generate_module_docs(self) -> None:
        """Generate documentation for all modules.

        :returns: None
        :rtype: None
        """
        meta_files: list[Path] = self.find_meta_files("modules")
        template: Template = self.jinja_env.get_template("module.md.j2")

        modules_output_dir: Path = self.output_dir / "docs" / "modules"
        # Remove pages of modules that were deleted since the last build
        shutil.rmtree(modules_output_dir, ignore_errors=True)
        modules_output_dir.mkdir(parents=True, exist_ok=True)

        module_index: list[dict[str, str]] = []

        for meta_file in meta_files:
            try:
                meta_data: dict[str, Any] = self.parse_meta_yml(meta_file)
                path_info: dict[str, str] = self.get_module_path_info(
                    meta_file, "modules"
                )

                # Process input/output for better template handling
                meta_data = self._process_meta_data(meta_data)

                # Create subdirectories as needed
                module_dir: Path = modules_output_dir / path_info["category"]
                if path_info["subcategory"]:
                    module_dir = module_dir / path_info["subcategory"]
                module_dir.mkdir(parents=True, exist_ok=True)

                # Get pipeline usage for this module
                pipelines = self.pipeline_usage_index.get(path_info["component"], [])

                # Generate documentation
                content: str = template.render(
                    meta=meta_data,
                    path_info=path_info,
                    module_type="module",
                    pipelines=pipelines,
                )

                output_file: Path = module_dir / f"{path_info['name']}.md"
                with open(output_file, "w") as f:
                    f.write(content)

                # Add to index
                module_index.append(
                    {
                        "name": meta_data.get("name", path_info["name"]),
                        "short_name": path_info["name"],
                        "description": meta_data.get("description", ""),
                        "path": str(output_file.relative_to(modules_output_dir)),
                        "category": path_info["category"],
                        "subcategory": path_info["subcategory"],
                    }
                )

                print(f"Generated docs for module: {path_info['full_path']}")

            except Exception as e:
                print(f"Error processing {meta_file}: {e}")

        # Generate modules index
        self._generate_index(module_index, modules_output_dir, "modules")

    def generate_subworkflow_docs(self) -> None:
        """Generate documentation for all subworkflows.

        :returns: None
        :rtype: None
        """
        meta_files: list[Path] = self.find_meta_files("subworkflows")
        template: Template = self.jinja_env.get_template("subworkflow.md.j2")

        subworkflows_output_dir: Path = self.output_dir / "docs" / "subworkflows"
        # Remove pages of subworkflows that were deleted since the last build
        shutil.rmtree(subworkflows_output_dir, ignore_errors=True)
        subworkflows_output_dir.mkdir(parents=True, exist_ok=True)

        subworkflow_index: list[dict[str, str]] = []

        for meta_file in meta_files:
            try:
                meta_data: dict[str, Any] = self.parse_meta_yml(meta_file)
                path_info: dict[str, str] = self.get_module_path_info(
                    meta_file, "subworkflows"
                )

                # Process input/output for better template handling
                meta_data = self._process_meta_data(meta_data)

                # Create subdirectories as needed
                subworkflow_dir: Path = subworkflows_output_dir / path_info["category"]
                if path_info["subcategory"]:
                    subworkflow_dir = subworkflow_dir / path_info["subcategory"]
                subworkflow_dir.mkdir(parents=True, exist_ok=True)

                # Get pipeline usage for this subworkflow
                pipelines = self.pipeline_usage_index.get(path_info["component"], [])

                # Generate documentation
                content: str = template.render(
                    meta=meta_data,
                    path_info=path_info,
                    module_type="subworkflow",
                    pipelines=pipelines,
                )

                output_file: Path = subworkflow_dir / f"{path_info['name']}.md"
                with open(output_file, "w") as f:
                    f.write(content)

                # Add to index
                subworkflow_index.append(
                    {
                        "name": meta_data.get("name", path_info["name"]),
                        "description": meta_data.get("description", ""),
                        "path": str(output_file.relative_to(subworkflows_output_dir)),
                        "category": path_info["category"],
                        "subcategory": path_info["subcategory"],
                    }
                )

                print(f"Generated docs for subworkflow: {path_info['full_path']}")

            except Exception as e:
                print(f"Error processing {meta_file}: {e}")

        # Generate subworkflows index
        self._generate_index(subworkflow_index, subworkflows_output_dir, "subworkflows")

    def _generate_index(
        self, items: list[dict[str, str]], output_dir: Path, item_type: str
    ) -> None:
        """Generate index page for modules or subworkflows.

        :param items: List of items to include in index
        :type items: list[dict[str, str]]
        :param output_dir: Output directory
        :type output_dir: Path
        :param item_type: 'modules' or 'subworkflows'
        :type item_type: str
        :returns: None
        :rtype: None
        """
        template: Template = self.jinja_env.get_template("index.md.j2")

        # Group by category
        categories: dict[str, list[dict[str, str]]] = {}
        for item in items:
            category: str = item["category"]
            if category not in categories:
                categories[category] = []
            categories[category].append(item)

        content: str = template.render(
            categories=categories, item_type=item_type, title=item_type.capitalize()
        )

        index_file: Path = output_dir / "index.md"
        with open(index_file, "w") as f:
            f.write(content)

    def generate_pipelines_docs(
        self,
        pipelines_config: Path,
        pipeline_usage_index: dict[str, list[dict[str, str]]],
    ) -> None:
        """Generate documentation page for pipelines.

        :param pipelines_config: Path to pipelines.yml configuration
        :type pipelines_config: Path
        :param pipeline_usage_index: Dictionary mapping module names to pipelines
        :type pipeline_usage_index: dict[str, list[dict[str, str]]]
        :returns: None
        :rtype: None
        """
        template: Template = self.jinja_env.get_template("pipelines.md.j2")

        # Load pipeline configuration
        with open(pipelines_config) as f:
            config = yaml.safe_load(f)
            pipelines_list = config.get("pipelines", [])

        # Build reverse index: pipeline -> modules/subworkflows
        pipeline_components: dict[str, dict[str, list[str]]] = {}
        for component_name, using_pipelines in pipeline_usage_index.items():
            for pipeline_info in using_pipelines:
                pipeline_name: str = pipeline_info["name"]
                if pipeline_name not in pipeline_components:
                    pipeline_components[pipeline_name] = {
                        "modules": [],
                        "subworkflows": [],
                    }

                # Determine if it's a module or subworkflow by checking path
                if (
                    self.modules_repo_path
                    / "modules"
                    / "ebi-metagenomics"
                    / component_name
                    / "meta.yml"
                ).exists():
                    pipeline_components[pipeline_name]["modules"].append(component_name)
                elif (
                    self.modules_repo_path
                    / "subworkflows"
                    / "ebi-metagenomics"
                    / component_name
                    / "meta.yml"
                ).exists():
                    pipeline_components[pipeline_name]["subworkflows"].append(
                        component_name
                    )

        # Enrich pipeline data with component information
        for pipeline in pipelines_list:
            pipeline_name = pipeline.get("name")
            if pipeline_name in pipeline_components:
                pipeline["modules"] = sorted(
                    pipeline_components[pipeline_name]["modules"]
                )
                pipeline["subworkflows"] = sorted(
                    pipeline_components[pipeline_name]["subworkflows"]
                )
            else:
                pipeline["modules"] = []
                pipeline["subworkflows"] = []

        # Pipelines are already sorted in pipelines.yml by:
        # 1. has_modules_json (true first)
        # 2. stars (descending)
        # So we keep the order from the YAML file

        # Generate documentation
        content: str = template.render(pipelines=pipelines_list)

        pipelines_output_dir: Path = self.output_dir / "docs"
        pipelines_output_dir.mkdir(parents=True, exist_ok=True)

        output_file: Path = pipelines_output_dir / "pipelines.md"
        with open(output_file, "w") as f:
            f.write(content)

        print(f"Generated pipelines documentation: {output_file}")


@click.command()
@click.option(
    "--modules-repo", required=True, help="Path to the nf-core modules repository"
)
@click.option(
    "--output-dir",
    default=".",
    help="Output directory for generated documentation (default: current directory)",
)
@click.option(
    "--pipelines-config",
    default="pipelines.yml",
    help="Path to pipelines configuration file (default: pipelines.yml)",
)
def main(modules_repo: str, output_dir: str, pipelines_config: str) -> None:
    """Generate MkDocs documentation from nf-core modules and subworkflows.

    :param modules_repo: Path to the nf-core modules repository
    :type modules_repo: str
    :param output_dir: Output directory for generated documentation
    :type output_dir: str
    :param pipelines_config: Path to pipelines configuration file
    :type pipelines_config: str
    :returns: None
    :rtype: None
    """
    # Build pipeline usage index
    print("Fetching pipeline usage information...")
    fetcher = PipelineUsageFetcher(pipelines_config)
    pipeline_usage_index = fetcher.build_module_usage_index()
    print(f"Found {len(pipeline_usage_index)} modules/subworkflows used by pipelines\n")

    parser: ModuleParser = ModuleParser(modules_repo, output_dir, pipeline_usage_index)

    print("Generating module documentation...")
    parser.generate_module_docs()

    print("Generating subworkflow documentation...")
    parser.generate_subworkflow_docs()

    print("Generating pipelines documentation...")
    parser.generate_pipelines_docs(Path(pipelines_config), pipeline_usage_index)

    print(f"Documentation generated in {output_dir}/docs/")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate MkDocs documentation from nf-core modules and subworkflows meta.yml files.
"""

from pathlib import Path
from typing import Any

import click
import yaml
from jinja2 import Environment, FileSystemLoader


class ModuleParser:
    """Parse nf-core module meta.yml files and generate documentation."""

    def __init__(self, modules_repo_path: str, output_dir: str):
        """Initialize the parser.

        :param modules_repo_path: Path to the nf-core modules repository
        :type modules_repo_path: str
        :param output_dir: Output directory for generated documentation
        :type output_dir: str
        """
        self.modules_repo_path = Path(modules_repo_path)
        self.output_dir = Path(output_dir)
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
        search_path = self.modules_repo_path / directory
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
        parts = meta_file.relative_to(self.modules_repo_path / base_type).parts[:-1]
        return {
            "full_path": "/".join(parts),
            "category": parts[0] if parts else "",
            "name": parts[-1] if parts else "",
            "subcategory": "/".join(parts[1:-1]) if len(parts) > 2 else "",
        }

    def _detect_component_types(self, components: list[str]) -> list[dict[str, str]]:
        """Detect whether components are modules or subworkflows.

        :param components: List of component names
        :type components: list[str]
        :returns: List of component info with type detection
        :rtype: list[dict[str, str]]
        """
        component_info = []

        for component in components:
            # Check if component exists as a subworkflow
            subworkflow_path = (
                self.modules_repo_path
                / "subworkflows"
                / "ebi-metagenomics"
                / component
                / "meta.yml"
            )

            if subworkflow_path.exists():
                component_type = "subworkflow"
            else:
                component_type = "module"

            component_info.append({"name": component, "type": component_type})

        return component_info

    def _process_meta_data(self, meta_data: dict[str, Any]) -> dict[str, Any]:
        """Process meta data to handle list structures in input/output.

        :param meta_data: Raw meta data from YAML
        :type meta_data: dict[str, Any]
        :returns: Processed meta data with flattened input/output structures
        :rtype: dict[str, Any]
        """
        processed = meta_data.copy()

        # Process tools section if it's a list
        if "tools" in processed and isinstance(processed["tools"], list):
            processed_tools = {}
            for item in processed["tools"]:
                if isinstance(item, dict):
                    processed_tools.update(item)
            processed["tools"] = processed_tools

        # Process input section if it's a list
        if "input" in processed and isinstance(processed["input"], list):
            processed_input = {}
            for item in processed["input"]:
                if isinstance(item, list):
                    # Handle nested list structure (like diamond/blastp)
                    for sub_item in item:
                        if isinstance(sub_item, dict):
                            processed_input.update(sub_item)
                elif isinstance(item, dict):
                    # Handle flat dict structure (like antismash)
                    processed_input.update(item)
            processed["input"] = processed_input

        # Process output section - handle both list and dict formats
        if "output" in processed:
            if isinstance(processed["output"], list):
                # Handle list format (like antismash)
                processed_output = {}
                for item in processed["output"]:
                    if isinstance(item, dict):
                        processed_output.update(item)
                processed["output"] = processed_output
            elif isinstance(processed["output"], dict):
                # Output is already a dict, but check for nested structures
                processed_output = {}
                for key, value in processed["output"].items():
                    if isinstance(value, list) and len(value) > 0:
                        # Handle list of dicts in output
                        if isinstance(value[0], dict):
                            # Take the first dict as the main definition
                            processed_output[key] = value[0]
                        else:
                            processed_output[key] = value
                    else:
                        processed_output[key] = value
                processed["output"] = processed_output

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
        meta_files = self.find_meta_files("modules")
        template = self.jinja_env.get_template("module.md.j2")

        modules_output_dir = self.output_dir / "docs" / "modules"
        modules_output_dir.mkdir(parents=True, exist_ok=True)

        module_index = []

        for meta_file in meta_files:
            try:
                meta_data = self.parse_meta_yml(meta_file)
                path_info = self.get_module_path_info(meta_file, "modules")

                # Process input/output for better template handling
                meta_data = self._process_meta_data(meta_data)

                # Create subdirectories as needed
                module_dir = modules_output_dir / path_info["category"]
                if path_info["subcategory"]:
                    module_dir = module_dir / path_info["subcategory"]
                module_dir.mkdir(parents=True, exist_ok=True)

                # Generate documentation
                content = template.render(
                    meta=meta_data, path_info=path_info, module_type="module"
                )

                output_file = module_dir / f"{path_info['name']}.md"
                with open(output_file, "w") as f:
                    f.write(content)

                # Add to index
                module_index.append(
                    {
                        "name": meta_data.get("name", path_info["name"]),
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
        meta_files = self.find_meta_files("subworkflows")
        template = self.jinja_env.get_template("subworkflow.md.j2")

        subworkflows_output_dir = self.output_dir / "docs" / "subworkflows"
        subworkflows_output_dir.mkdir(parents=True, exist_ok=True)

        subworkflow_index = []

        for meta_file in meta_files:
            try:
                meta_data = self.parse_meta_yml(meta_file)
                path_info = self.get_module_path_info(meta_file, "subworkflows")

                # Process input/output for better template handling
                meta_data = self._process_meta_data(meta_data)

                # Create subdirectories as needed
                subworkflow_dir = subworkflows_output_dir / path_info["category"]
                if path_info["subcategory"]:
                    subworkflow_dir = subworkflow_dir / path_info["subcategory"]
                subworkflow_dir.mkdir(parents=True, exist_ok=True)

                # Generate documentation
                content = template.render(
                    meta=meta_data, path_info=path_info, module_type="subworkflow"
                )

                output_file = subworkflow_dir / f"{path_info['name']}.md"
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
        self, items: list[dict], output_dir: Path, item_type: str
    ) -> None:
        """Generate index page for modules or subworkflows.

        :param items: List of items to include in index
        :type items: list[dict]
        :param output_dir: Output directory
        :type output_dir: Path
        :param item_type: 'modules' or 'subworkflows'
        :type item_type: str
        :returns: None
        :rtype: None
        """
        template = self.jinja_env.get_template("index.md.j2")

        # Group by category
        categories = {}
        for item in items:
            category = item["category"]
            if category not in categories:
                categories[category] = []
            categories[category].append(item)

        content = template.render(
            categories=categories, item_type=item_type, title=item_type.capitalize()
        )

        index_file = output_dir / "index.md"
        with open(index_file, "w") as f:
            f.write(content)


@click.command()
@click.option(
    "--modules-repo", required=True, help="Path to the nf-core modules repository"
)
@click.option(
    "--output-dir",
    default=".",
    help="Output directory for generated documentation (default: current directory)",
)
def main(modules_repo: str, output_dir: str):
    """Generate MkDocs documentation from nf-core modules and subworkflows.

    :param modules_repo: Path to the nf-core modules repository
    :type modules_repo: str
    :param output_dir: Output directory for generated documentation
    :type output_dir: str
    :returns: None
    :rtype: None
    """
    parser = ModuleParser(modules_repo, output_dir)

    print("Generating module documentation...")
    parser.generate_module_docs()

    print("Generating subworkflow documentation...")
    parser.generate_subworkflow_docs()

    print(f"Documentation generated in {output_dir}/docs/")


if __name__ == "__main__":
    main()

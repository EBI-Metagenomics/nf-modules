# Contributing to Documentation

This guide explains how to contribute to the EBI Metagenomics nf-core modules documentation.

## Overview

The documentation is automatically generated from `meta.yml` files in the modules and subworkflows directories using MkDocs and custom Python scripts. The website is automatically deployed to GitHub Pages when changes are merged to the main branch.

## Quick Start

### Prerequisites

- [Task](https://taskfile.dev/installation/) - Task runner
- [uv](https://github.com/astral-sh/uv#installation) - Fast Python package installer

### Local Development

1. **Build documentation locally:**

   ```bash
   task docs-build
   ```

2. **Serve documentation locally:**

   ```bash
   task docs-serve
   ```

   This will start a local development server at http://127.0.0.1:8000

3. **Clean documentation:**
   ```bash
   task docs-clean
   ```

## Documentation Structure

```
docs/
├── docs/                    # Generated documentation content
│   ├── index.md            # Homepage
│   ├── modules/            # Module documentation (auto-generated)
│   └── subworkflows/       # Subworkflow documentation (auto-generated)
├── templates/              # Jinja2 templates for documentation generation
│   ├── index.md.j2        # Index page template
│   ├── module.md.j2       # Module page template
│   └── subworkflow.md.j2  # Subworkflow page template
├── generate_docs.py        # Documentation generation script
├── mkdocs.yml             # MkDocs configuration
└── pyproject.toml         # Python dependencies
```

## Contributing

### Module/Subworkflow Documentation

Documentation is automatically generated from `meta.yml` files. To improve documentation:

1. **Edit meta.yml files** in the respective module/subworkflow directories
2. **Key sections that affect documentation:**
   - `name`: Module/subworkflow name
   - `description`: Brief description
   - `keywords`: Tags for categorization
   - `tools`: Tool information with homepage and documentation links
   - `input`: Input specifications
   - `output`: Output specifications
   - `authors`: Author information
   - `maintainers`: Maintainer information

### Template Customization

To modify how documentation is generated:

1. **Edit templates** in `docs/templates/`:

   - `module.md.j2`: Template for individual module pages
   - `subworkflow.md.j2`: Template for individual subworkflow pages
   - `index.md.j2`: Template for index pages

2. **Template variables available:**
   - `meta`: All data from meta.yml
   - `path_info`: Path information (category, name, etc.)
   - `module_type`: "module" or "subworkflow"

### Website Configuration

To modify the website appearance or behavior:

1. **Edit `mkdocs.yml`** for:

   - Site metadata
   - Theme configuration
   - Navigation structure
   - Plugins and extensions

2. **Custom CSS** can be added to `docs/docs/stylesheets/extra.css`

### Generation Script

The `generate_docs.py` script handles:

  - Parsing meta.yml files
  - Processing input/output structures
  - Rendering Jinja2 templates
  - Creating directory structures

To modify generation logic, edit this script.

## Workflow

### Automatic Deployment

  - Documentation is automatically rebuilt and deployed when:
  - Changes are pushed to the main branch
  - meta.yml files are modified
  - Documentation files are changed
  - Workflow file is updated

### Manual Deployment

For testing or manual deployment:

```bash
# Generate documentation
cd docs
uv run python generate_docs.py --modules-repo .. --output-dir .

# Build and serve locally
uv run mkdocs serve

# Build for production
uv run mkdocs build
```

## Best Practices

### meta.yml Guidelines

1. **Use clear, concise descriptions**
2. **Include comprehensive input/output specifications**
3. **Provide tool homepage and documentation links**
4. **Keep author and maintainer information up to date**
5. **Use consistent naming conventions**

### Template Guidelines

1. **Keep templates readable and maintainable**
2. **Use consistent formatting across templates**
3. **Handle missing or optional data gracefully**
4. **Follow Jinja2 best practices**

## Troubleshooting

### Common Issues

1. **Documentation not updating:**

   - Check that meta.yml files are valid YAML
   - Verify GitHub Actions workflow completed successfully
   - Clear browser cache

2. **Local development issues:**

   - Ensure uv and dependencies are installed
   - Run `task docs-clean` and rebuild
   - Check Python virtual environment

3. **Template errors:**
   - Validate Jinja2 syntax
   - Check that all referenced variables exist
   - Test with minimal data

### Getting Help

- Check existing issues in the repository
- Review the MkDocs documentation: https://www.mkdocs.org/
- Consult the Material theme docs: https://squidfunk.github.io/mkdocs-material/

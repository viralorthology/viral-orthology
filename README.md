# Viral Orthology

A modular bioinformatics pipeline for automated orthology inference across viral genomes using complementary sequence, profile, synteny, compositional, and structural evidence.

> [!WARNING]
> Viral Orthology is currently undergoing a major rewrite focused on improving
> code quality, modularity, maintainability, and test coverage. The rewrite is
> still under active development and **is not yet ready for installation or use**.

## Rewrite Progress

- [x] Core objects
- [x] Core engines
- [ ] Core pipeline
- [ ] Enrichment modules
- [ ] Analysis modules

## Key Features

- Automated orthology inference across viral genomes
- Iterative refinement of orthologous groups
- Support for large viral genome datasets
- Configurable parameters for integrated bioinformatics tools
- Independent enrichment and analysis modules
- Synteny detection based on gene-order conservation
- Amino acid composition analysis
- Protein secondary and tertiary structure analysis

## Workflow

ViralOrthology consists of a core orthology inference pipeline followed by optional enrichment and analysis modules.
The core pipeline processes the input sequences, generates initial orthologous groups, and iteratively refines them.
The resulting orthologous groups can then be processed by enrichment modules, which use complementary evidence such as synteny conservation, amino acid composition, and protein structural information.
Finally, analysis modules can be used to examine the resulting groups, evaluate their quality, and extract additional information from the inferred relationships.

![Viral Orthology pipeline](workflow.png)

## Installation

### Requirements

- Linux
- Conda (or Miniconda)

### Install

```bash
curl -fsSL https://raw.githubusercontent.com/viralorthology/viral-orthology/refs/heads/main/setup.sh | bash
```

The required dependencies are installed automatically through the Conda environment.

### Verify installation

```bash
viralorthology -h
```

> [!NOTE]
> The commands below describe the planned command-line interface and may not be functional until the rewrite is complete.

## Quick Start

### 1. Prepare the input

Create a file named `ids.txt` containing the GenBank accession IDs of the viral genomes to be analyzed:

```text
NC_XXXXX
NC_XXXXX
NC_XXXXX
```

### 2. Retrieve sequences

```bash
viralorthology -download_seqs
```

This command retrieves the required genomic and protein sequence data and prepares the input for downstream analysis.

### 3. Run the core pipeline

```bash
viralorthology -pipeline
```

The core workflow performs the initial orthology inference and iterative refinement of the resulting orthologous groups.

### 4. Run enrichment and analysis modules

Once the orthologous groups have been generated, independent enrichment and analysis modules can be executed as needed.

For example:

```bash
viralorthology -synteny
```

```bash
viralorthology -composition ...
```

## Documentation

Detailed documentation is available on [viralorthology.github.io](https://viralorthology.github.io).

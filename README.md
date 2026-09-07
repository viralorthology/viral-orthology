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

## Key Features

* Automated orthology inference across viral genomes
* Iterative refinement of orthologous groups
* Support for large viral genome datasets
* Configurable parameters for integrated bioinformatics tools
* Independent enrichment and assessment modules
* Synteny detection based on gene-order conservation
* Amino acid composition analysis
* Protein secondary and tertiary structure analysis

## Architecture

The pipeline is organized into two major components:

* **Core pipeline** — generates and iteratively refines orthologous groups.
* **Enrichment modules** — independently analyze the resulting groups using synteny, amino acid composition, and protein structure.

Enrichment modules are designed to operate independently from the core orthology inference workflow.

## Workflow

The analysis begins with the automated core workflow, which processes the input sequences and generates orthologous groups. Independent downstream modules can then be applied to enrich or assess the resulting groups.

## Installation

### Requirements

* Linux
* Conda (or Miniconda)

### Install

```bash
# Installation instructions will be added here
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

### 4. Run downstream modules

Once the orthologous groups have been generated, independent downstream modules can be executed according to the desired analysis.

For example:

```bash
viralorthology -synteny
```

```bash
viralorthology -composition ...
```

## Documentation

Detailed documentation is available on [viralorthology.github.io](https://viralorthology.github.io).

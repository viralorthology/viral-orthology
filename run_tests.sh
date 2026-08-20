#!/usr/bin/env bash
set -e

ruff format --check src
ruff check src
mypy --strict src
PYTHONPATH=src pytest --cov=src

#!/usr/bin/env bash
set -e

PYTHONPATH=src pytest
ruff check src
mypy --strict src

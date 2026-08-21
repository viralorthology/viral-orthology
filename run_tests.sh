#!/usr/bin/env bash
set -e

ruff format --check src
ruff check src
mypy --strict --disable-error-code=no-untyped-call src
PYTHONPATH=src pytest --cov=src

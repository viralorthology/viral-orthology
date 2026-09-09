#!/bin/bash
set -euo pipefail

local_bin_path="$HOME/.local/bin"
viralorthology_files_path="$HOME/.viralorthology"
mkdir -p "$local_bin_path"
rm -rf "$viralorthology_files_path"
mkdir -p "$viralorthology_files_path"

# check for conda
if ! command -v conda >/dev/null 2>&1; then
    echo "Error: conda is required but was not found." >&2
    exit 1
fi

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

# download viralorthology
wget https://github.com/viralorthology/viral-orthology/archive/refs/heads/main.zip -O "$tmp_dir/viral-orthology.zip"
unzip -q "$tmp_dir/viral-orthology.zip" -d "$tmp_dir"

# install 3rd party libs
# cp -rp libs/* "$local_bin_path" # TODO

# install viralorthology
cp -rp "$tmp_dir/viral-orthology-main/src/"* "$viralorthology_files_path"

# make or update conda env
if conda run -n viralorthology true 2>/dev/null; then
    echo "Conda environment 'viralorthology' already exists. Skipping."
else
    conda env create -n viralorthology -f "$tmp_dir/viral-orthology-main/environment.yml"
fi

# make launcher
cat << "EOF" > "$local_bin_path/viralorthology"
#!/usr/bin/env bash

conda run --no-capture-output -n viralorthology python "$HOME/.viralorthology/viralorthology.py" "$@"
EOF
chmod +x "$local_bin_path/viralorthology"

# download ORFfinder
if [[ ! -f "$local_bin_path/ORFfinder" ]]; then
    wget -O "$local_bin_path/ORFfinder.gz" https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
    gzip -d "$local_bin_path/ORFfinder.gz"
    chmod +x "$local_bin_path/ORFfinder"
fi

# download InterProScan
if [[ ! -f "$HOME/.interproscan/interproscan.sh" ]]; then
    rm -rf "$HOME/.interproscan"
    wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.78-109.0/interproscan-5.78-109.0-64-bit.tar.gz -O "$tmp_dir/interproscan.tar.gz"
    tar -xzf "$tmp_dir/interproscan.tar.gz" -C "$tmp_dir"
    mv "$tmp_dir/interproscan-5.78-109.0" "$HOME/.interproscan"
    ln -s "$HOME/.interproscan/interproscan.sh" "$local_bin_path/interproscan"
fi

# check if home/.local/bin is in path
if [[ ":$PATH:" != *":$local_bin_path:"* ]]; then
    echo "Warning: $local_bin_path is not in your PATH."
    echo "Add the following line to your shell configuration:"
    echo "export PATH=\"$local_bin_path:\$PATH\""
fi

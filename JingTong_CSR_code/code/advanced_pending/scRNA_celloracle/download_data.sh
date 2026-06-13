#!/bin/bash
# OPTIONAL — only if you did NOT upload the GSE160756/*.loom data with the package.
# Fetches GSE160756 scRNA (7 .loom: 3 NP + 2 CEP + 2 AF) directly on the server.
set -e
mkdir -p GSE160756 && cd GSE160756
echo "downloading GSE160756_RAW.tar (~433 MB) ..."
wget -c "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE160nnn/GSE160756/suppl/GSE160756_RAW.tar"
tar -xf GSE160756_RAW.tar
gunzip -f *.loom.gz
rm -f GSE160756_RAW.tar
echo "done:"; ls -la *.loom

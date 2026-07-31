#!/bin/bash

set -e

ACCESSION="SRR23609077"

mkdir -p input

prefetch "$ACCESSION"

fasterq-dump "$ACCESSION" \
    --split-files \
    --threads 4 \
    --outdir input

gzip input/${ACCESSION}_1.fastq
gzip input/${ACCESSION}_2.fastq


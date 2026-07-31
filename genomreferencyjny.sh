#!/bin/bash

mkdir -p genom

wget -O genom/NC_045512.2.fasta \
"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_045512.2&rettype=fasta&retmode=text"


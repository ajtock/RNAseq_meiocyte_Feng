#!/bin/bash

source activate RNAseq_mapping
snakemake -p --cores 16
conda deactivate

#!/usr/bin/env bash


snakemake --cores 16 --use-conda --conda-frontend mamba -p --rerun-incomplete --cache


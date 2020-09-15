#!/bin/bash

python3 generation.py
bwa index ref.fa
bwa mem ref.fa reads/reads_substitution.fasta -o sam_files/subs.sam
bwa mem ref.fa reads/reads_all_events.fasta -o sam_files/all_events.sam

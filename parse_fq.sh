#!/bin/bash

paste -d $'\t' <(paste -d $'\t' - - - - < sample_1.fastq) <(paste -d $'\t' - - - - < sample_2.fastq) | awk -v FS=$'\t' '$6~/^CCTTGGACACCCGAGAATTCCA/' | cut -f 1-4 -d $'\t' | tr $'\t' '\n' > sample_1.fq.filtered

paste -d $'\t' <(paste -d $'\t' - - - - < sample_1.fastq) <(paste -d $'\t' - - - - < sample_2.fastq) | awk -v FS=$'\t' '$6~/^CCTTGGACACCCGAGAATTCCA/' | cut -f 5-8 -d $'\t' | tr $'\t' '\n' > sample_2.fq.filtered

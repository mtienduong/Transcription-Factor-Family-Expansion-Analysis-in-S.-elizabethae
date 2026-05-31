## Overview

This project is a follow-up analysis to our transcription factor (TF) family enrichment study in Synalpheus elizabethae and Synalpheus brooksi.

In the previous analysis, transcripts were assigned to transcription factor families using the CrusTF database, and TF family distributions were compared between the full transcriptome and differentially expressed genes (DEGs) identified between queens and workers.

Several TF families appeared to contain multiple transcripts. However, an increase in transcript number can arise from two different biological processes:

Alternative splicing, where a single gene produces multiple transcript isoforms.
Gene family expansion, where multiple related genes are present in the genome.

Because gene family expansion and contraction have been associated with the evolution of eusociality, this analysis aims to determine whether highly represented TF clusters are primarily composed of transcript isoforms or distinct gene copies.

## Biological Question

When a transcription factor family contains many transcripts, are those transcripts:

Multiple isoforms of the same Trinity gene?
Multiple related genes belonging to the same TF family?
A combination of both?

Answering this question helps distinguish transcript-level diversity from potential gene family expansion.

## Project Summary
This analysis:

Clusters TF-associated transcripts using MMseqs2 at 50% sequence identity.

Identifies clusters containing large numbers of transcripts.

Counts the number of transcripts present in each cluster.

Estimates the number of unique Trinity genes represented within each cluster.
Determines whether cluster size is driven by isoform diversity or multiple genes.

## Tools Used

MMseqs2: Sequence clustering

Apptainer/Singularity: Containerized MMseqs2 execution

AWK: Cluster parsing and summary statistics

## Workflow
### 1. Extract TF-associated transcripts
TF-associated transcripts were identified in the previous CrusTF analysis.
Ex: S.brooksi. For DEG-specific analyses, transcript sequences were extracted using:
```bash
seqtk subseq \
S.brooksi_transcriptome_assembly.fasta \
S.brooksi_TF_DEGs.txt \
> S.brooksi_TF_DEGs_sequences.fasta
```
### 2. Cluster sequences with MMseqs2
Sequences were clustered at 50% amino acid identity.
```bash
mmseqs easy-cluster \
INPUT_FASTA \
OUTPUT_NAME \
TEMP_DIR \
--min-seq-id 0.5 \
--threads 36
```
This produces a cluster file where:
ClusterRepresentative     ClusterMember

Each row links a transcript to its assigned cluster.

### 3. Identify large clusters

Count the number of transcripts associated with each cluster representative:
```bash
cut -f1 S.brooksi_all_transcript_clusters_50.tsv \
| sort \
| uniq -c \
| awk '{print $2 "\t" $1}'
```

### 4. Estimate Trinity isoform abundance
To determine the expected range of transcript isoforms per Trinity gene:
```bash
grep ">" Trinity.fasta \
| cut -d' ' -f1 \
| sed 's/_i[0-9]*//' \
| sort \
| uniq -c \
| sort -nr \
| head
```
Example output:

42 >TRINITY_DN6645_c1_g1

37 >TRINITY_DN369_c0_g1

33 >TRINITY_DN2509_c3_g1

### 5. Distinguish Isoforms from Gene Family Expansion

Transcript IDs were collapsed to the Trinity gene level by removing isoform suffixes (_i#).

```bash
awk '{
rep=$1 # representative transcript (cluster rep)
member=$2 # member transcript

# Extract gene part by removing _i# from transcript IDs
rep_gene=rep
member_gene=member
sub(/_i[0-9]+$/,"",rep_gene)
sub(/_i[0-9]+$/,"",member_gene)

# Count total members
total[rep]++

# Count isoforms (same gene as rep) vs different gene
if(member_gene==rep_gene) iso[rep]++
else dup[rep]++
}
END
```

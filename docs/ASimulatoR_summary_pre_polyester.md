# ASimulatoR Summary - Pre-Polyester Steps

**Date:** 2026-02-11
**Focus:** Simplified summary of ASimulatoR pipeline steps before polyester simulation

---

## Overview

This document provides a simplified explanation of the main steps in ASimulatoR up to the point where polyester is called. The pipeline transforms a genome annotation (GTF/GFF3) into transcript variants with alternative splicing events, preparing everything needed for RNA-seq read simulation.

---

## Step 1: Exon Superset Creation

**Goal**: Create a "template" containing all possible exons for each gene

**Process**:
1. Collect all exons from all transcripts of a gene
2. Merge overlapping exons together (using interval reduction)
3. Fix the order based on strand orientation (+ or -)
4. Calculate transcript positions (1, 2, 3, ...) in sequence order

**Key Formula**:
```
S_i = reduce(union of all exons for gene g_i)
```

**Example**: If a gene has three annotated transcripts with exons at positions [1-100, 201-300], [50-150, 201-300], and [1-100, 151-250], the exon superset would be [1-150, 151-300] after merging all overlapping regions.

---

## Step 2: Assign Events to Genes

**Goal**: Decide which genes get which alternative splicing events

**Process**:

1. **Gene Length Distribution**: Count how many exons each gene has
2. **Cumulative Distribution**: Calculate how many genes have at least x exons
3. **Two Assignment Modes**:
   - **Probability mode**: Each compatible gene has a random chance to get each event
   - **Frequency mode**: Exact percentage of genes get each event (deterministic)
4. **Resource Constraint Check**: Ensure enough genes have enough exons for each event type
5. **Iterative Adjustment**: If constraints violated, reduce total genes and retry

**Key Formulas**:
```
Genes per event (frequency mode): n_c = floor(p_c × N)

Minimum exons needed for combination: r_c = sum(r_e for all events) - (number_of_events - 1)

Genes with enough exons: C_x = count of genes with ≥ x exons
```

---

## Step 3: Map Events to Exons

**Goal**: Place events on specific exons within each selected gene

**Process**:
1. Randomly pick an event from the set to place
2. Randomly pick a starting position in the available exon list
3. Assign that event to a contiguous block of exons
4. Split remaining exons into left and right segments
5. Partition remaining events randomly between segments
6. Recursively repeat for each segment
7. If constraints violated, backtrack and retry (up to 10 attempts)
8. Fall back to greedy placement if random placement fails

**Key Constraints**:
```
Non-overlap: No exon can be used by more than one event
Minimum size: Each event type needs enough exons (es needs 3, mes needs 4+, etc.)
Feasibility: Events must fit within available exon "real estate"
```

---

## Step 4: Build Variant Transcripts

**Goal**: Create the actual spliced transcript by applying event rules

**Process**:
- Mark which exons to include/exclude based on event type
- **Exon skipping (es)**: Skip middle exon in the assigned block
- **Multiple exon skipping (mes)**: Skip intermediate exons in the assigned block
- **Alternative first exon (afe)**: Randomly include one of first two exons
- **Alternative last exon (ale)**: Randomly include one of last two exons
- **Mutually exclusive exons (mee)**: Randomly include one of two middle exons
- **Intron retention (ir)**: Merge two exons, keep the intron between them
- **Alternative splice sites (a3, a5)**: Modify exon boundaries within interior

**Transcript coordinate calculation**:
```
tr_start[1] = 1
tr_start[j] = tr_end[j-1] + 1  (for j > 1)
tr_end[j] = tr_start[j] + exon_length - 1
```

---

## Step 5: Create Event Annotations

**Goal**: Record what happened for downstream evaluation

**Process**: For each event, save:
- Event type (es, mes, ir, a3, a5, mee, afe, ale)
- Variant transcript ID
- Template transcript ID
- Genomic coordinates (chromosome positions)
- Transcriptomic coordinates (position within transcript)

This annotation serves as ground truth for evaluating AS detection tools.

---

## Step 6: Generate Fold Changes Matrix

**Goal**: Set up differential expression between experimental groups

**Is this step part of polyester?**

**No, this is NOT part of polyester.** This step is performed by **ASimulatoR** before calling `polyester::simulate_experiment()`. ASimulatoR generates the fold changes matrix and passes it as a parameter to polyester. Polyester then uses this matrix to adjust read counts per group.

---

### When No Differential Expression is Set

If differential expression is NOT explicitly set, the algorithm **automatically** decides whether to introduce differential expression:

**The Decision Process**:

For each gene with multiple transcripts:

```
Random number between 0 and 1 ──► If < 0.5 (50% chance):
                                  ✓ No differential expression
                                  (All fold changes = 1)
                               Else (50% chance):
                                  ✓ Add isoform switch
                                  (One transcript gets 2x change)
```

**Key insight**: Even if you don't request differential expression, there's a 50% chance genes with multiple transcripts will have an isoform switch. This is designed to simulate realistic biological variation.

---

### Detailed Algorithm Breakdown

#### Determine Number of Groups
```
nr_groups = 2                      (if num_reps is NULL, default)
nr_groups = length(num_reps)       (if num_reps is specified)
```

#### For Each Gene

**Case A: Gene has only 1 transcript**
```
Result matrix:
  Group1  Group2  ...
    1       1      1
```
**Explanation**: With only one transcript, differential expression within the gene is impossible.

**Case B: Gene has multiple transcripts**

**Path 1** (50% probability): No differential expression
```
Example for gene with 3 transcripts, 2 groups:
  Group1  Group2
    1       1      (transcript 1)
    1       1      (transcript 2)
    1       1      (transcript 3)
```

**Path 2** (50% probability): Isoform switch
```
Step 1: Create first column (Group1)
  Transcript 1: 2  (selected randomly for 2x change)
  Transcript 2: 1
  Transcript 3: 1

Step 2: For additional groups (Group2, etc.):
  - Shuffle the first column randomly
  - Ensure it's not identical to the first column
  - This means the isoform switch occurs in different groups

Example outcome:
  Group1  Group2
    2       1      (transcript 1: high in group1, low in group2)
    1       2      (transcript 2: low in group1, high in group2)
    1       1      (transcript 3: same in both groups)
```

---

### What the Fold Changes Matrix Means

The matrix tells polyester how many reads to generate for each transcript in each group:

```
fc_t = fold change for transcript t in current group

Interpretation:
- fc_t = 1:  Baseline expression
- fc_t = 2:  2x higher expression (double the reads)
- fc_t = 0.5: 2x lower expression (half the reads)
```

---

### The Formula in Context

```
Expected reads per transcript:
r_t = D × (l_t / total_length) × fc_t
```

**Breaking it down**:

| Component | Meaning |
|-----------|---------|
| `D` | Total sequencing depth (reads per sample) |
| `l_t / total_length` | Proportion of total transcriptome length represented by transcript t |
| `fc_t` | **Fold change multiplier from ASimulatoR's matrix** |

**Example calculation**:
```
D = 10,000,000 reads
l_t (transcript 1) = 1,000 bp
total_length = 100,000,000 bp (sum of all transcripts)
fc_t = 2 (isoform switch upregulated)

r_t = 10,000,000 × (1,000/100,000,000) × 2
    = 10,000,000 × 0.01 × 2
    = 200,000 reads
```

Without the fold change (fc_t = 1), this transcript would get 100,000 reads. With fc_t = 2, it gets 200,000 reads.

---

## Summary Table

| Step | Goal | Key Output |
|------|------|------------|
| 1: Exon Superset Creation | Merge all exons per gene | Exon supersets with transcript coordinates |
| 2: Assign Events to Genes | Decide which genes get events | Gene-to-event assignments |
| 3: Map Events to Exons | Place events on specific exons | Event-to-exon mappings |
| 4: Build Variant Transcripts | Create actual spliced transcripts | Variant transcripts (GTF) |
| 5: Create Event Annotations | Record ground truth | Event annotation table |
| 6: Generate Fold Changes Matrix | Set up differential expression | Fold changes matrix (passed to polyester) |

---

## Most Important Formulas Summary

1. **Exon superset**:
   ```
   S_i = reduce(union of exons)
   ```

2. **Minimum exons for events**:
   ```
   r_c = sum(r_e) - (|E| - 1)
   ```
   (Subtract (|E|-1) because events can share boundary exons)

3. **Transcript coordinates**:
   ```
   tr_end[j] = tr_start[j] + width - 1
   ```

4. **Read expectation**:
   ```
   r_t = D × (l_t / total_length) × fc_t
   ```

5. **Resource check**:
   ```
   If C_{r_c} < n_c → need to reduce N (total genes)
   ```
   Where C_{r_c} = genes with at least r_c exons, n_c = genes allocated to event c

---

## Key Parameters

| Parameter | Meaning |
|-----------|---------|
| `event_probs` | Probabilities or frequencies for each AS event type |
| `probs_as_freq` | If TRUE, use deterministic frequency mode |
| `max_genes` | Maximum genes to include (NULL = all available) |
| `seq_depth` | Total reads per sample (D) |
| `num_reps` | Number of replicates per group |
| `fold_changes` | Matrix generated by ASimulatoR, passed to polyester |

---

## Output Files (Pre-Polyester)

| File | Description |
|------|-------------|
| `splicing_variants.gtf` | All transcript variants (templates + AS variants) |
| `event_annotation.tsv` | Ground truth for AS events |
| `splicing_variants_novel.gtf` | Subset for novel event detection |
| Exon-junction table | Coverage information for validation |

These files, along with the fold changes matrix, are then passed to polyester for read simulation.

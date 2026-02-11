# ASimulatoR Mathematical Representation

**Date:** 2026-02-05  
**Focus:** Bioinformatic science details of the ASimulatoR pipeline up to transcript construction and count matrix creation

## 1. Pipeline Overview

ASimulatoR is a splice-aware RNA-Seq data simulator that generates gold standard datasets for evaluating alternative splicing (AS) analysis tools. The pipeline consists of three main stages:

1. **Exon Superset Creation**: Merges all annotated exons from genome annotation into unspliced templates
   - Input: GTF/GFF3 genome annotation file
   - Output: Exon supersets (RDA format) with transcriptomic coordinates
   - This creates a "template" for each gene containing all possible exons

2. **Transcript Variant Creation**: Generates transcript variants with user-defined AS event distributions
   - Input: Exon supersets, event probabilities/frequencies
   - Output: GTF/GFF3 annotation files, event annotation table, exon-junction coverage table
   - This stage constructs variants by applying AS events to templates

3. **RNA-Seq Read Simulation**: Simulates RNA-Seq reads using a modified version of polyester
   - Input: GTF annotation, chromosome FASTA files
   - Output: FASTQ files, count matrix, transcript information
   - This stage generates sequencing reads from the transcript variants

**Data Flow**:
```
GTF/GFF3 → Exon Supersets → Transcript Variants + Event Annotation → Polyester → FASTQ + Counts
```

**Key Parameters**:
- `event_probs`: Probabilities or relative frequencies for each AS event type
- `probs_as_freq`: If TRUE, treat event_probs as relative frequencies (deterministic allocation)
- `max_genes`: Maximum number of genes to include (NULL = all available)
- `seq_depth`: Total reads per sample
- `num_reps`: Number of replicates per group

**Output Files**:
- `splicing_variants.gtf`: All transcript variants (templates + AS variants)
- `event_annotation.tsv`: Ground truth for AS events
- `splicing_variants_novel.gtf`: Subset annotation for novel event detection
- FASTQ files: Simulated sequencing reads
- `counts.txt`: Read counts per transcript per sample

## 2. Mathematical Notation

### 2.1 Basic Definitions

- $G = \{g_1, g_2, \ldots, g_n\}$: Set of genes in the genome annotation
- For each gene $g_i$:
  - $T_i = \{t_{i1}, t_{i2}, \ldots, t_{ik_i}\}$: Set of transcript variants
  - $E_i = \{e_{i1}, e_{i2}, \ldots, e_{im_i}\}$: Set of exons from all transcript variants (superset)
  - Each exon $e_{ij} = (chr, start_{ij}, end_{ij}, strand)$: Genomic coordinates
  - $m_i$: Number of exons in exon superset of gene $g_i$

**Coordinate Systems**:
- **Genomic coordinates**: Absolute positions on the chromosome (1-based)
- **Transcriptomic coordinates**: Relative positions within the mature transcript (1-based, 5'→3')

### 2.2 AS Event Types

ASimulatoR supports eight AS event types with minimum exon requirements:

Let $\mathcal{E} = \{es, mes, ir, a3, a5, mee, afe, ale\}$ be the set of AS event types with minimum exon requirements:
- $r_{es} = 3$: Exon Skipping (skipped exon requires flanking exons)
- $r_{mes} = k + 2$ where $k \geq 2$: Multiple Exon Skipping  
- $r_{ir} = 2$: Intron Retention (retained intron between two exons)
- $r_{a3} = 2$: Alternative 3' Splice Site (alternative acceptor site)
- $r_{a5} = 2$: Alternative 5' Splice Site (alternative donor site)
- $r_{mee} = 4$: Mutually Exclusive Exons (choose one of two middle exons)
- $r_{afe} = 2$: Alternative First Exon (swap first two exons)
- $r_{ale} = 2$: Alternative Last Exon (swap last two exons)

## 3. Exon Superset Creation

### 3.1 Exon Merging Algorithm

For each gene $g_i$, create an exon superset $S_i$ by:

1. **Exon Collection**: Gather all exons from all transcript variants:
   $$E_i = \bigcup_{t \in T_i} E_{it}$$
   where $T_i$ is the set of transcript variants for gene $g_i$, and $E_{it}$ are the exons of transcript $t$.

2. **Exon Reduction**: Merge overlapping exons using GenomicRanges::reduce():
   $$S_i = \text{reduce}(E_i)$$
   This creates a set of non-overlapping genomic intervals representing all possible exon positions. The `reduce()` function merges overlapping or adjacent intervals, creating a minimal set of non-overlapping intervals that cover all original exons.

3. **Strand Orientation**: For negative strand genes, reverse the order:
   $$S_i = \begin{cases}
   S_i & \text{if strand} = + \\
   \text{rev}(S_i) & \text{if strand} = -
   \end{cases}$$
   This ensures consistent 5'→3' ordering regardless of genomic strand.

4. **Transcriptomic Coordinates**: Compute transcriptomic coordinates:
   $$\text{tr\_start}_{ij} = \begin{cases}
   1 & \text{if } j = 1 \\
   \text{tr\_end}_{i(j-1)} + 1 & \text{otherwise}
   \end{cases}$$
   $$\text{tr\_end}_{ij} = \text{tr\_start}_{ij} + \text{width}(e_{ij}) - 1$$
   where $\text{width}(e_{ij}) = \text{end}_{ij} - \text{start}_{ij} + 1$.

### 3.2 Mathematical Representation

The exon superset for gene $g_i$ is:
$$S_i = \{(e_{ij}, \text{tr\_start}_{ij}, \text{tr\_end}_{ij}, \text{gene\_exon\_number}_j) \mid j = 1,\ldots,|S_i|\}$$
where $|S_i|$ is the number of exons in the superset.

**Key Insight**: The exon superset represents the *union* of all possible exonic regions for a gene, creating a "template" from which alternative splicing variants can be constructed by including or excluding specific exons.

## 4. Transcript Variant Creation

### 4.1 Event Probability Model

User specifies event probabilities $P = \{p_e \mid e \in \mathcal{E}'\}$ where $\mathcal{E}' \subseteq \mathcal{E}$ can include event combinations (e.g., $ES,IR$).

Two modes:
1. **Probability mode**: $p_e$ is probability for each compatible exon superset
2. **Frequency mode**: $p_e$ is relative frequency, requiring $\sum_{e \in \mathcal{E}'} p_e \leq 1$

### 4.2 Gene Selection Algorithm (Relative Frequency Mode, No max_genes)

**When `probs_as_freq = TRUE` and `max_genes = NULL`**, the algorithm uses all available compatible exon supersets with deterministic event assignment based on relative frequencies. This algorithm ensures that exactly the specified relative frequencies of AS events are maintained across the generated dataset.

#### 4.2.1 Initial Setup

1. **Gene Length Distribution**: For all genes $G$, compute:
   $$L = \{l_i = |S_i| \mid g_i \in G\}$$
   where $l_i$ is the number of exons in gene $g_i$'s superset. Single-exon genes cannot host AS events since at least two exons are required for any splicing event.

2. **Cumulative Distribution**: Create cumulative counts of genes with at least $x$ exons:
   $$C_x = \sum_{i=1}^n \mathbb{I}(l_i \geq x)$$
   This is computed as `cumsum(rev(table(factor(gene_lengths, levels = 1:max(gene_lengths)))))` - a reverse cumulative sum of the frequency table. The value $C_x$ represents the number of genes that have $x$ or more exons, which is the pool from which events requiring $x$ exons can be assigned.

3. **Available Genes**: Set $N = \sum_{i=1}^n \mathbb{I}(l_i > 1)$ (genes with >1 exon, since single-exon genes cannot have AS events). This is the initial total number of genes to be used in the simulation.

#### 4.2.2 Deterministic Event Assignment (Relative Frequency Mode)

For relative frequencies $p_c$ (where $\sum_c p_c \leq 1$):

1. **Allocate Genes per Event**: For each event combination $c$:
   $$n_c = \lfloor p_c \cdot N \rfloor$$
   Genes are deterministically assigned to events based on these counts. For example, if $p_{es} = 0.1$ and $N = 1000$, then exactly 100 genes will have exon skipping events. This deterministic allocation ensures reproducibility and precise control over event distributions.

2. **Create Assignment Matrix**: Build binary matrix $A$ of size $N \times |\mathcal{E}'|$ where:
   - First $n_{c_1}$ genes get event $c_1$
   - Next $n_{c_2}$ genes get event $c_2$
   - etc.
   - Remaining genes get no event (if $\sum_c p_c < 1$)

    In code: `construct_all <- sapply(names(event_probs), function(event) { rep(names(event_probs) == event, floor(event_probs[[event]] * nr_genes)) })`

3. **Flatten to List**: Convert matrix to list of event assignments per gene. Each gene gets exactly one event combination (or no event) when using relative frequency mode.

**Key Difference from Probability Mode**: In relative frequency mode, the total number of genes assigned to each event is fixed by the specified frequencies, ensuring that the final dataset has exactly the desired proportions of AS events. In probability mode, each compatible gene has an independent probability of hosting each event type, introducing randomness in the final event counts.

#### 4.2.3 Resource Constraint Checking

For each event combination $c$, compute minimum exon requirement:
$$r_c = \sum_{e \in c} r_e - (|c| - 1)$$

The formula subtracts $(|c| - 1)$ because when multiple events are assigned to the same gene, they can share boundary exons. For example, if a gene has events {es, a3} where es requires 3 exons and a3 requires 2 exons, the total minimum exons needed is 3 + 2 - 1 = 4, not 5, because the a3 event's first exon can be the same as the es event's third exon.

Check if sufficient genes meet this requirement:
$$\text{If } C_{r_c} < n_c \text{ for any } c$$

Where $C_{r_c}$ is the number of genes with at least $r_c$ exons (obtained from the cumulative distribution).

**Example**: If 100 genes are allocated to exon skipping ($r_{es} = 3$) but only 80 genes have 3 or more exons, the constraint is violated and the algorithm must reduce the total number of genes.

#### 4.2.4 Iterative Adjustment

If insufficient genes are available for the requested event distribution, the algorithm iteratively reduces the total number of genes until constraints can be satisfied:

1. Compute missing genes: $m_c = n_c - C_{r_c}$
   - $m_c > 0$ indicates how many more genes with at least $r_c$ exons are needed
   - $m_c \leq 0$ means there are enough genes available for event $c$

2. Reduce total genes proportionally:
   $$N \leftarrow N - \lceil \frac{N}{n_c} \cdot m_c \rceil$$
   This reduction is proportional to ensure that the relative frequencies of events are maintained as closely as possible while satisfying all constraints.

3. Recompute allocations and check again

**Loop continues until all event requirements can be satisfied with available gene resources.** The iteration stops when for all event combinations $c$, we have $C_{r_c} \geq n_c$.

#### 4.2.5 Gene Selection

After event assignments are finalized:
1. Sort genes by their event requirements (highest first)
2. For each gene assignment, randomly select a compatible gene:
   $$\text{gene}_i \sim \text{Uniform}(\{g \in G \mid l_g \geq r_{c_i} \text{ and not yet selected}\})$$
3. Mark selected genes as unavailable for subsequent selections

### 4.3 Event Assignment to Exons

For each selected gene $g_i$ with assigned events $E_i \subseteq \mathcal{E}'$:

#### 4.3.1 Event Exon Mapping Function

Define $f: E_i \to 2^{\{1,\ldots,|S_i|\}}$ mapping events to exon indices, with constraints:
- $|f(es)| = 3$: middle exon skipped, flanking exons included
- $|f(mes)| = k + 2$ where $k \geq 2$: $k$ middle exons skipped
- $|f(ir)| = 2$: intron between these exons retained
- $|f(a3)| = 2$: alternative acceptor site in second exon
- $|f(a5)| = 2$: alternative donor site in first exon
- $|f(mee)| = 4$: mutually exclusive middle exons (positions 2 and 3)
- $|f(afe)| = 2$: first two exons swapped
- $|f(ale)| = 2$: last two exons swapped

#### 4.3.2 Recursive Assignment Algorithm

The `.assign_events()` function implements a **divide-and-conquer recursive algorithm** that partitions exon indices among competing AS events. This algorithm ensures that multiple events can coexist on the same transcript without overlapping in their exon usage.

**Algorithm Design Philosophy**: The algorithm uses a stochastic approach with a fallback to greedy assignment. It attempts to place events randomly to create diverse transcript structures. If random placement fails after `max_attempts` tries, it falls back to placing the next event at the first available position, ensuring that the algorithm always finds a solution if one exists.

**Input**:
- $V = \{v_1, v_2, \ldots, v_m\}$: Available exon indices (1-based positions in template)
- $R = \{(e_1, r_1), (e_2, r_2), \ldots, (e_k, r_k)\}$: Events with their minimum exon requirements $r_e$

**Algorithm**:

1. **Base Case**: If $R = \emptyset$, return $\emptyset$
   - No events remain to assign, so return an empty assignment

2. **Recursive Step** (attempts up to `max_attempts = 10` times):
    a. **Random Event Selection**: Choose event $e \sim \text{Uniform}(R)$
       - Random selection promotes diversity in event placement
       - Events compete for exon space in a non-deterministic manner

    b. **Random Position Selection**: If attempt ≤ max_attempts:
       $$s \sim \text{Uniform}\{1, 2, \ldots, |V| - r_e + 1\}$$
       Otherwise (greedy fallback): $s = 1$
       - The greedy fallback ensures termination when random placement consistently fails
       - Position $s$ is the starting exon index for the current event

    c. **Assign Exons**: $f(e) = \{v_s, v_{s+1}, \ldots, v_{s+r_e-1}\}$
       - The event is assigned to a contiguous block of $r_e$ exons
       - This reflects the biological reality that AS events involve contiguous exon regions

    d. **Split Remaining Exons**:
       $$V_1 = \{v_1, \ldots, v_{s-1}\}$$
       $$V_2 = \{v_{s+r_e}, \ldots, v_{|V|}\}$$
       - The exons used by event $e$ are removed from the pool
       - Remaining exons are split into left and right segments for recursive assignment

    e. **Partition Remaining Events**: Randomly split $R \setminus \{e\}$ between $V_1$ and $V_2$:
       $$\text{split} \sim \text{Binomial}(|R|-1, 0.5)$$
       Events assigned to $V_1$ and $V_2$ must satisfy:
       $$\sum_{e \in R_1} r_e - (|R_1| - 1) \leq |V_1|$$
       $$\sum_{e \in R_2} r_e - (|R_2| - 1) \leq |V_2|$$
       - The constraint checks ensure that each partition has enough exons for its assigned events
       - If the random partition violates constraints, the algorithm retries from step 2b
       - This backtracking mechanism ensures only feasible assignments are accepted

    f. **Recursive Calls**:
       $$f_1 = \text{assign\_events}(V_1, R_1)$$
       $$f_2 = \text{assign\_events}(V_2, R_2)$$

    g. **Success Condition**: If $|f_1| + |f_2| = |R| - 1$, return $\{f(e)\} \cup f_1 \cup f_2$
       - This condition verifies that all remaining events were successfully assigned

**Key Constraints**:
1. **Non-overlap**: $\forall e_i, e_j \in R, i \neq j: f(e_i) \cap f(e_j) = \emptyset$
2. **Compact assignment**: Events assigned to contiguous exon blocks
3. **Resource feasibility**: Event combinations must fit within available exons

**Example**: For events `{es:3, a3:2}` on 5-exon gene:
- Possible assignment: `es → {1,2,3}`, `a3 → {4,5}`
- Not allowed: `es → {1,2,3}`, `a3 → {2,3}` (overlap)
- Not allowed: `es → {1,2,3}`, `a3 → {5}` (insufficient exons for a3)

**Mathematical Insight**: This is a **constrained combinatorial allocation problem** where events compete for exon "real estate". The algorithm ensures biological plausibility by:
1. Maintaining exon contiguity within events
2. Preventing physically impossible overlapping events
3. Respecting minimum exon requirements for each event type

### 4.4 Transcript Construction

For template $S_i$ and event assignments $f$, construct variant transcript $V_i$:

#### 4.4.1 Exon Selection Function

Define inclusion vector $I \in \{0,1\}^{|S_i|}$ where $I_j = 1$ if exon $j$ included in variant.

Application rules for each event type:

The `apply_exon_events()` function transforms inclusion vector $I$:

1. **Exon Skipping (es)**: Skip middle exon
   $$I_{f(es)_2} \leftarrow 0$$
   - Exon skipping requires three exons: flanking exons (positions 1 and 3) are retained, and the middle exon (position 2) is skipped
   - This simulates the most common AS event type where a cassette exon is excluded from the mature transcript

2. **Multiple Exon Skipping (mes)**: Skip intermediate exons
   $$I_{f(mes)_k} \leftarrow 0 \quad \forall k \in \{2,\ldots,|f(mes)|-1\}$$
   - The number of skipped exons is randomly determined: $k \sim \text{Uniform}\{2, \ldots, \text{available\_exons} + 2 - \min\_nr\_exons\}$
   - The first and last exons in the assigned block are retained as flanking exons
   - This allows simulation of more complex skipping events involving multiple consecutive exons

3. **Alternative First Exon (afe)**: Swap first two exons
   $$I_1, I_2 \leftarrow I_2, I_1$$
   - This simulates alternative promoter usage where the first exon of the transcript is replaced
   - Only one of the first two exons is included in the final variant (the other is excluded)
   - The exon to include is randomly selected

4. **Alternative Last Exon (ale)**: Swap last two exons
   $$I_{m-1}, I_m \leftarrow I_m, I_{m-1}$$
   - This simulates alternative polyadenylation site usage
   - Only one of the last two exons is included in the final variant

5. **Mutually Exclusive Exons (mee)**: Choose one of two middle exons
   $$I_{f(mee)_2}, I_{f(mee)_3} \leftarrow I_{f(mee)_3}, I_{f(mee)_2}$$
   - This simulates mutually exclusive exon selection where exactly one of two adjacent exons is included
   - The choice is random, and the order swap ensures exactly one exon is retained

6. **Intron Retention (ir)**: Merge two exons (handled separately)
   - The intron between two exons is retained as part of the mature transcript
   - This is handled by merging the two flanking exons into a single exon with extended boundaries

7. **Alternative Splice Sites**: Modify boundaries (handled separately):
   - For $a3$ (Alternative 3' Splice Site): $end_{f(a3)_2} \leftarrow new\_end \in (start_{f(a3)_2} + 1, end_{f(a3)_2} - 1)$
     - The acceptor site of the second exon is shifted, shortening or lengthening the exon
     - The new position is randomly sampled from within the exon's interior
   - For $a5$ (Alternative 5' Splice Site): $start_{f(a5)_1} \leftarrow new\_start \in (start_{f(a5)_1} + 1, end_{f(a5)_1} - 1)$
     - The donor site of the first exon is shifted
     - If both a5 and a3 events are present on the same exon, the new positions are constrained to avoid overlap

**Coordinate Transformation for Merged Exons (IR)**:
For intron retention, the two exons are merged:
$$\text{new\_start} = \min(\text{start}_{exon1}, \text{start}_{exon2})$$
$$\text{new\_end} = \max(\text{end}_{exon1}, \text{end}_{exon2})$$
The retained intron is:
$$(\text{end}_{exon1} + 1, \text{start}_{exon2} - 1)$$

**Transcriptomic Coordinate Recalculation**:
After applying events, the transcriptomic coordinates must be recalculated:
$$\text{tr\_start}_{ij} = \begin{cases}
1 & \text{if } j = 1 \\
\text{tr\_end}_{i(j-1)} + 1 & \text{otherwise}
\end{cases}$$
$$\text{tr\_end}_{ij} = \text{tr\_start}_{ij} + \text{width}(e_{ij}) - 1$$
This ensures that the transcriptomic coordinates reflect the actual sequence of the variant transcript.

#### 4.4.2 Variant Transcript

The variant transcript is:
$$V_i = \{S_{ij} \mid I_j = 1\} \quad \text{with modified boundaries for A3SS/A5SS/IR}$$

### 4.5 Event Annotation Generation

The `get_event_annotation()` function creates records for each event, capturing both genomic and transcriptomic coordinates. This annotation serves as the ground truth for evaluating AS detection tools.

**Annotation Structure**:
Each annotation record contains:
- `event_annotation`: Event type (e.g., "es", "mes", "ir")
- `variant`: Transcript ID of the variant transcript
- `template`: Transcript ID of the reference (template) transcript
- `genomic_start`, `genomic_end`: Chromosome coordinates
- `transcriptomic_start`, `transcriptomic_end`: Positions within the transcript sequence

#### 4.5.1 Detailed Event Annotation Rules

1. **Exon Skipping (es)**:
   - Genomic: $(\text{start}_{skipped}, \text{end}_{skipped})$ - coordinates of the skipped exon in the template
   - Transcriptomic: $(\text{tr\_start}_{skipped}, \text{tr\_end}_{skipped})$ - positions in the template transcript
   - Note: The skipped exon is present in the template but absent in the variant

2. **Multiple Exon Skipping (mes)**:
   - Genomic: $(\text{start}_{skipped_1},\text{end}_{skipped_1}),\ldots,(\text{start}_{skipped_k},\text{end}_{skipped_k})$
   - Transcriptomic: $(\text{tr\_start}_{skipped_1},\text{tr\_end}_{skipped_1}),\ldots$
   - Multiple exons are skipped; their coordinates are comma-separated

3. **Intron Retention (ir)**:
   - Genomic: $(\text{end}_{upstream} + 1, \text{start}_{downstream} - 1)$
     - The intron spans from the end of the upstream exon to the start of the downstream exon
   - Transcriptomic: $(tr\_start_{retained}, tr\_end_{retained})$ where:
     $$tr\_start_{retained} = tr\_start_{upstream} + \Delta$$
     $$\Delta = \begin{cases}
     \text{end}_{upstream} - \text{end}_{ri} & \text{if strand} = - \\
     \text{start}_{ri} - \text{start}_{upstream} & \text{otherwise}
     \end{cases}$$
   - The $\Delta$ term accounts for strand-specific orientation when calculating transcriptomic coordinates

4. **Alternative Splice Sites**:
   - **A3SS (Alternative 3' Splice Site)**:
     - Genomic: $C_g = (\min(new\_end, old\_end) + 1, \max(new\_end, old\_end))$
       - This captures the region between the old and new acceptor sites
     - Transcriptomic: $C_t = (tr\_start_{exon}, tr\_start_{exon} + |new\_end - old\_end| - 1)$
       - The affected region starts at the exon's transcriptomic start and spans the distance between splice sites
   - **A5SS (Alternative 5' Splice Site)**:
     - Genomic: $C_g = (\min(new\_start, old\_start), \max(new\_start, old\_start) - 1)$
     - Transcriptomic: $C_t = (tr\_end_{exon} - |new\_start - old\_start| + 1, tr\_end_{exon})$
       - For negative strand genes, the order of coordinates is reversed

5. **Mutually Exclusive Exons (mee)**: Two annotation entries are created:
   - Entry 1: Coordinates of the included exon in the variant
   - Entry 2: Coordinates of the excluded exon in the template
   - This allows direct comparison between the two mutually exclusive choices

6. **Alternative First Exon (afe)**: Two annotation entries are created:
   - Entry 1: Coordinates of the variant's first exon
   - Entry 2: Coordinates of the template's first exon
   - This captures the alternative promoter usage

7. **Alternative Last Exon (ale)**: Two annotation entries are created:
   - Entry 1: Coordinates of the variant's last exon
   - Entry 2: Coordinates of the template's last exon
   - This captures the alternative polyadenylation site usage

**Strand Considerations**:
For negative strand genes, the genomic coordinate calculations are inverted because transcription occurs in the 3'→5' direction on the reference genome. The transcriptomic coordinates, however, are always 5'→3' relative to the mature transcript.

## 5. Preparation for Polyester Simulation

### 5.1 Transcript Count Calculation

Before calling polyester, ASimulatoR calculates the number of transcripts per gene from the exon-junction table:

$$tr\_per\_gene = \{t_g \mid g \in G\}$$
where $t_g = |\{v \mid v \in V, v.gene\_id = g\}|$ is the number of transcript variants for gene $g$.

This count includes both the template transcript (no AS events) and all generated variant transcripts.

### 5.2 Fold Changes Matrix Generation

For evaluating differential splicing tools, ASimulatoR introduces fold changes between groups. The fold changes matrix $F$ is generated using a stochastic algorithm that creates isoform switches between groups.

**Definitions**:
- $n_{groups}$: Number of experimental groups (default: 2)
- $n_{reps}$: Vector specifying number of replicates per group
- $F$: Fold changes matrix of size $(\sum_{g \in G} t_g) \times n_{groups}$

**Algorithm**:

1. **Determine number of groups**:
   $$n_{groups} = \begin{cases}
   2 & \text{if } num\_reps \text{ is } NULL \\
   length(num\_reps) & \text{otherwise}
   \end{cases}$$

2. **For each gene** with $t_g$ transcript variants:
   - If $t_g = 1$ (single transcript): Create a fold changes matrix with no differential expression
     $$F_g = \begin{pmatrix} 1 & 1 & \ldots & 1 \end{pmatrix}$$
     - All transcripts have fold change of 1 (no change between groups)

   - If $t_g > 1$ (multiple transcripts): With probability 0.5:
     - **No differential splicing**: All transcripts have fold change of 1
       $$F_g = \mathbf{1}_{t_g \times n_{groups}}$$

     - **Isoform switch**: Introduce a 2-fold change for one randomly selected transcript
       $$F_g[:, 1] = \begin{pmatrix} 2 \\ 1 \\ \vdots \\ 1 \end{pmatrix}$$
       - The first transcript gets a 2-fold increase in group 1
       - All other transcripts remain at baseline (fold change = 1)

     - For additional groups ($n_{groups} > 1$): Create column vectors by sampling the first column
       $$F_g[:, i] = \text{sample}(F_g[:, 1]) \quad \text{for } i = 2, \ldots, n_{groups}$$
       - Constraint: Ensure $F_g[:, i] \neq F_g[:, 1]$ for at least one transcript (avoid identical groups)
       - This constraint is enforced by resampling if the vectors are identical

3. **Combine matrices**:
   $$F = \begin{pmatrix} F_{g_1} \\ F_{g_2} \\ \vdots \\ F_{g_{|G|}} \end{pmatrix}$$

**Interpretation**:
- Fold change of 2 means the transcript is twice as abundant in that group compared to the reference
- Fold change of 1 means no differential expression
- The matrix is row-indexed by transcript ID, following the order in which polyester will process the transcripts

**Example**: For a gene with 3 transcripts across 2 groups:
$$F_g = \begin{pmatrix} 1 & 1 \\ 1 & 2 \\ 1 & 1 \end{pmatrix}$$
This indicates that transcript 2 is differentially expressed, with a 2-fold increase in group 2 compared to group 1.

**Current Limitation**: As noted in the code, "Currently, ASimulatoR introduces random isoform switches. Those can be retraced in the sim_tx_info.txt file written by polyester. We plan on improving this in the future."

### 5.3 Polyester Parameters Configuration

ASimulatoR sets specific polyester parameters before simulation:

| Parameter | Value | Description |
|-----------|-------|-------------|
| `gtf` | `outdir/splicing_variants.gtf` | GTF file containing all transcripts |
| `seqpath` | `input_dir` | Directory containing chromosome FASTA files |
| `outdir` | `outdir` | Output directory for FASTQ files |
| `meanmodel` | `TRUE` | Reads per transcript as function of transcript length |
| `frag_GC_bias` | `NULL` | GC bias disabled (not supported in ASimulatoR) |
| `fold_changes` | $F$ | Fold changes matrix for differential splicing |
| `strand_specific` | `TRUE` | Strand-specific simulation |
| `exon_junction_table` | Exon-junction table | For coverage calculation |
| `exon_junction_coverage` | `TRUE` | Enable coverage output |

**Transcript Ordering**: Polyester processes transcripts in the order they appear in the GTF file. The fold changes matrix must be ordered consistently with this transcript order.

## 6. RNA-Seq Read Simulation (Overview)

### 6.1 Polyester Integration

ASimulatoR uses a modified version of polyester with parameters:

- **Sequencing depth**: $D$ reads per sample
- **Read length**: $L$ bases
- **Error rate**: $\epsilon$
- **Fragment length distribution**: $FL \sim \text{Normal}(\mu_f, \sigma_f^2)$
- **Strand-specific**: Boolean

**Modifications to Polyester**:
ASimulatoR extends polyester with two technical bias simulations:

1. **PCR Duplicates**:
   - User specifies fraction $d$ of reads to duplicate
   - Number of duplications per read follows $k \sim \text{Poisson}(\lambda=1)$
   - Default $\lambda = 1$ ensures most reads have 0 or 1 duplication
   - Duplicated reads are marked in the output

2. **Adapter Contamination**:
   - For fragments shorter than read length $L$, append adapter sequence
   - This simulates sequencing beyond the fragment end where adapter is read
   - Adapter sequence is user-configurable

### 6.2 Read Generation Model

For each transcript $t$ with length $l_t$:

1. **Expected reads calculation** (with meanmodel enabled):
   $$r_t = D \cdot \frac{l_t}{\sum_{t'} l_{t'}} \cdot fc_t$$
   where:
   - $D$ is the sequencing depth (total reads per sample)
   - $l_t$ is the length of transcript $t$
   - $\sum_{t'} l_{t'}$ is the total length of all transcripts
   - $fc_t$ is the fold change for transcript $t$ in the current group

   This model ensures that longer transcripts receive proportionally more reads, reflecting the biological reality that longer transcripts generate more fragments during library preparation.

2. **Fragment generation**: Sample $r_t$ fragment start positions:
   $$s_i \sim \text{Uniform}(1, l_t - FL + 1) \quad \text{for } i = 1, \ldots, r_t$$
   - Each fragment has length $FL \sim \text{Normal}(\mu_f, \sigma_f^2)$
   - Fragment start positions are uniformly distributed across the transcript

3. **Read simulation**: For each fragment, generate paired-end reads with:
   - **Read 1**: Sequence from $s_i$ to $s_i + L - 1$
   - **Read 2**: Sequence from $s_i + FL - L$ to $s_i + FL - 1$
   - **Sequencing errors**: Each base error with probability $\epsilon$
   - **Quality scores**: Phred-scaled based on error probability:
     $$Q = -10 \cdot \log_{10}(\epsilon)$$

4. **Strand-specificity**: When enabled, reads are oriented correctly:
   - Read 1 maps to forward strand (5'→3' direction)
   - Read 2 maps to reverse strand (3'→5' direction on reference)
   - This preserves strand information for splice junction detection

### 6.3 Technical Biases

1. **PCR Duplicates**: Fraction $d$ of reads duplicated, duplication count $\sim \text{Poisson}(\lambda=1)$
   - PCR amplification during library preparation creates duplicate copies of some fragments
   - This bias affects quantitative analyses of read counts
   - Marking duplicates allows downstream tools to account for this bias

2. **Adapter Contamination**: For fragments shorter than $L$, append adapter sequence
   - When fragment length < read length, sequencing reads into the adapter
   - This occurs more frequently with degraded RNA or small RNA species
   - Trimming is required during read processing to remove adapter sequences

## 7. Illustrative Examples

### 7.1 Example: Exon Skipping Event

Consider a gene with the following exon superset (genomic coordinates):

```
Exon 1: chr1:1000-1200 (+)
Exon 2: chr1:3000-3200 (+)  <-- will be skipped
Exon 3: chr1:5000-5200 (+)
```

**Template transcript coordinates**:
```
Genomic:  1000-1200 | 3000-3200 | 5000-5200
Tr-coord: 1-201      | 202-402   | 403-603
```

**Exon skipping (es) event**:
- Event assigned to exons {1, 2, 3}
- Middle exon (exon 2) is skipped

**Variant transcript**:
```
Genomic:  1000-1200 | 5000-5200
Tr-coord: 1-201      | 202-402
```

**Event annotation**:
```
event:    es
variant:  gene1_es
template: gene1_template
genomic:  3000-3200 (skipped exon)
tr-coord: 202-402 (position in template)
```

### 7.2 Example: Multiple Event Combination

Consider a gene with 5 exons hosting both exon skipping (es) and alternative 3' splice site (a3):

**Exon superset**:
```
Exon 1: chr1:1000-1200
Exon 2: chr1:3000-3200  <-- skipped in ES
Exon 3: chr1:5000-5300  <-- a3 affects end position
Exon 4: chr1:7000-7200
Exon 5: chr1:9000-9200
```

**Event assignment**:
- es → {1, 2, 3} (skips exon 2)
- a3 → {3, 4} (modifies end of exon 3)

**Minimum exon calculation**:
- es requires 3 exons: $r_{es} = 3$
- a3 requires 2 exons: $r_{a3} = 2$
- Combined: $\min\_exons(\{es, a3\}) = 3 + 2 - 1 = 4$
  - The -1 accounts for shared boundary (exon 3 is part of both events)

**Possible recursive assignment**:
```
Iteration 1:
  - Randomly select event: es
  - Randomly select position: 1
  - Assign es → {1, 2, 3}
  - Remaining exons: V1 = {}, V2 = {4, 5}
  - Split remaining events: a3 → V2 (only option since V1 empty)
  - Assign a3 → {4, 5}
Result: es → {1,2,3}, a3 → {4,5}
```

**Variant construction**:
1. Apply es: skip exon 2 → include {1, 3, 4, 5}
2. Apply a3: modify end of exon 3 (in variant, this is exon 2)
   - Original end: 5300
   - New end: 5100 (randomly selected within exon)
3. Final variant:
```
Exon 1: 1000-1200 (tr-coord: 1-201)
Exon 2: 5000-5100 (tr-coord: 202-301) [modified end]
Exon 3: 7000-7200 (tr-coord: 302-402)
Exon 4: 9000-9200 (tr-coord: 403-503)
```

### 7.3 Example: Fold Changes Matrix

Consider 3 genes with different numbers of transcripts across 2 groups:

**Gene 1**: 2 transcripts (t1, t2)
- With 50% probability: isoform switch
- Fold changes: t1 = 2, t2 = 1
- Matrix: `[[2, 1], [1, 1]]` or `[[1, 1], [1, 2]]` (random)

**Gene 2**: 1 transcript (t3)
- No differential expression possible
- Matrix: `[[1, 1]]`

**Gene 3**: 3 transcripts (t4, t5, t6)
- With 50% probability: no differential expression
- Matrix: `[[1, 1], [1, 1], [1, 1]]`
- Or isoform switch: `[[2, 1], [1, 2], [1, 1]]` (example)

**Combined fold changes matrix**:
```
       Group1  Group2
t1        2       1
t2        1       1
t3        1       1
t4        1       1
t5        1       1
t6        1       1
```

This matrix tells polyester to generate 2x as many reads for transcript t1 in Group1 compared to Group2.

### 7.4 Example: Resource Constraint Adjustment

Assume:
- 1000 genes available
- 200 genes have ≥ 3 exons
- Requested frequencies: es = 0.3 (300 genes), mes = 0.2 (200 genes)

**Initial allocation**:
- es: $n_{es} = 0.3 \times 1000 = 300$ genes
- mes: $n_{mes} = 0.2 \times 1000 = 200$ genes

**Constraint check**:
- es requires 3 exons: need 300 genes with ≥ 3 exons
- mes requires 4 exons: need 200 genes with ≥ 4 exons
- Available: 200 genes with ≥ 3 exons, only 150 with ≥ 4 exons

**Violation detected**: mes needs 200 genes but only 150 available

**Adjustment**:
- Missing genes for mes: $m_{mes} = 200 - 150 = 50$
- Reduce total: $N \leftarrow 1000 - \lceil (1000/200) \times 50 \rceil = 750$

**Recalculate**:
- es: $n_{es} = 0.3 \times 750 = 225$ genes
- mes: $n_{mes} = 0.2 \times 750 = 150$ genes

**Verify**: Now 150 genes available for mes (exactly what's needed)

The algorithm iterates until all constraints are satisfied, ensuring the final dataset is feasible.

## 8. Output Files and Data Structures

### 8.1 GTF/GFF3 Export Files

ASimulatoR generates several annotation files that are used as inputs to polyester and for downstream analysis:

**1. splicing_variants.gtf**: Primary annotation file
- Contains all transcript variants including templates and AS variants
- Features include: `gene`, `transcript`, `exon`, `junction`, `ri` (retained intron)
- Junction features represent exon-exon boundaries for splice junction detection
- Retained intron features indicate introns that are included in the mature transcript

**2. splicing_variants.gff3**: GFF3 format equivalent
- Optional output (controlled by `write_gff` parameter)
- Contains `ID` and `Parent` attributes for hierarchical relationships
- Useful for tools that require GFF3 format

**3. splicing_variants_novel.gtf**: Subset annotation
- Optional output (controlled by `novel_variants` parameter)
- Contains only a subset (1 - novel_variants) of all variants
- Used to evaluate AS tool performance on "novel" (unannotated) events
- Novel events are randomly selected from the full variant set

### 8.2 Event Annotation Table

The `event_annotation.tsv` file contains ground truth information about all AS events:

**Columns**:
- `event_annotation`: Event type (es, mes, ir, a3, a5, mee, afe, ale)
- `variant`: Transcript ID of the variant containing the event
- `template`: Transcript ID of the reference (template) transcript
- `genomic_start`, `genomic_end`: Chromosome coordinates
- `transcriptomic_start`, `transcriptomic_end`: Positions within transcript

**Uses**:
- Ground truth for evaluating AS detection tool accuracy
- Precision and recall calculations
- Validation of event boundary predictions

### 8.3 Exon-Junction Coverage Table

When `exon_junction_coverage = TRUE`, ASimulatoR returns a data table with:
- Exon coordinates and coverage information
- Junction coordinates
- Retained intron coordinates

This table is used internally for coverage calculation and validation.

### 8.4 Polyester Output Files

Polyester generates the following output files:

**1. FASTQ files**:
- `sample_01_read1.fastq`, `sample_01_read2.fastq`: Paired-end reads for sample 1
- Additional samples follow the same naming pattern
- FASTQ files contain read sequences and Phred quality scores

**2. Count matrix**:
- `counts.txt`: Read counts for each transcript in each sample
- Rows: Transcripts (ordered by appearance in GTF file)
- Columns: Samples
- This matrix can be used for differential expression analysis

**3. Transcript information**:
- `sim_tx_info.txt`: Metadata about simulated transcripts
- Contains fold changes, expression levels, and other simulation parameters
- Allows tracing of which transcripts were differentially expressed

## 9. Mathematical Summary

### 9.1 Key Equations

1. **Exon superset size**: $|S_i| = |\text{reduce}(\bigcup_{t \in T_i} E_{it})|$

2. **Minimum exons for event combination**:
   $$\min\_exons(E) = \sum_{e \in E} r_e - (|E| - 1)$$
   where $r_e$ is minimum exon requirement for event $e$
   - The $(|E| - 1)$ term accounts for shared boundary exons when multiple events co-occur

3. **Event assignment probability**: For probability mode:
   $$P(\text{gene } g_i \text{ has event } e) = p_e \cdot \mathbb{I}(|S_i| \geq r_e)$$
   - Each compatible gene independently has probability $p_e$ of hosting event $e$

4. **Variant construction**:
   $$V_i = \text{apply\_events}(S_i, f, I)$$
   where $f$ is event-exon mapping and $I$ is inclusion vector

5. **Transcriptomic coordinate calculation**:
   $$\text{tr\_start}_{ij} = \begin{cases}
   1 & \text{if } j = 1 \\
   \text{tr\_end}_{i(j-1)} + 1 & \text{otherwise}
   \end{cases}$$
   $$\text{tr\_end}_{ij} = \text{tr\_start}_{ij} + \text{width}(e_{ij}) - 1$$

6. **Read count expectation (with meanmodel)**:
   $$r_t = D \cdot \frac{l_t}{\sum_{t'} l_{t'}} \cdot fc_t$$

### 9.2 Algorithm Complexity

- **Exon superset creation**: $O(\sum_i |T_i| \cdot \log |E_i|)$
  - Dominated by interval reduction operation on each gene's exon set

- **Event assignment**: $O(N \cdot |\mathcal{E}'| \cdot \text{max\_attempts})$
  - $N$: Number of genes
  - $|\mathcal{E}'|$: Number of event types
  - max_attempts: Maximum number of placement attempts per event (default: 10)

- **Variant construction**: $O(\sum_i |S_i| \cdot |E_i|)$
  - Linear in total number of exons across all genes

- **Fold changes generation**: $O(|G| \cdot \text{avg\_transcripts\_per\_gene} \cdot n_{groups})$
  - Depends on number of genes, average transcripts per gene, and number of groups

## 10. Biological Interpretation and Use Cases

### 10.1 Biological Accuracy

The mathematical model captures several key biological aspects of alternative splicing:

1. **Exon combinatorics**: All possible exon combinations from annotated transcripts
   - The exon superset represents the union of all possible exonic regions for a gene
   - This allows reconstruction of any annotated transcript variant

2. **AS event constraints**: Biological feasibility of event types based on exon count
   - Minimum exon requirements ensure that simulated events are biologically plausible
   - For example, exon skipping requires at least 3 exons (flanking exons + skipped exon)

3. **Strand awareness**: Proper orientation of negative strand genes
   - Coordinate transformations account for strand-specific transcription direction
   - Transcriptomic coordinates are always 5'→3' relative to the mature transcript

4. **Coordinate systems**: Both genomic and transcriptomic coordinates for AS analysis
   - Genomic coordinates enable mapping reads back to the reference genome
   - Transcriptomic coordinates facilitate splice junction annotation and analysis

### 10.2 Use Cases for AS Tool Evaluation

ASimulatoR's mathematical model enables several types of evaluation:

1. **Event type-specific performance**:
   - By controlling event distributions, users can evaluate tool performance on specific AS event types
   - Example: Create a dataset with 50% exon skipping to test ES detection sensitivity

2. **Sequencing depth analysis**:
   - Varying $D$ (sequencing depth) allows studying the effect of depth on detection sensitivity
   - Can determine minimum depth required for reliable AS event detection

3. **Technical bias assessment**:
   - PCR duplication and adapter contamination parameters allow studying tool robustness to technical artifacts
   - Tools can be evaluated on their ability to handle biased data

4. **Differential splicing evaluation**:
   - Fold changes matrix creation enables testing differential splicing detection tools
   - Ground truth fold changes allow calculation of sensitivity and specificity

5. **Novel event detection**:
   - The `novel_variants` parameter creates datasets where some events are not annotated
   - This allows evaluating tool performance on discovering novel AS events

### 10.3 Limitations and Future Improvements

1. **Fold changes generation**: Current implementation uses random isoform switches
   - Future improvements may allow user-specified fold change patterns
   - More sophisticated models of differential splicing could be incorporated

2. **Expression level control**: Current meanmodel is length-based only
   - Biological variation in transcript expression is not modeled
   - Future versions may incorporate tissue-specific expression patterns

3. **Event co-occurrence constraints**: Current model allows limited multi-event combinations
   - Biological constraints on which events can co-occur could be incorporated
   - More realistic models of complex multi-exon transcripts could be developed

### 10.4 Mathematical Validity

The mathematical model ensures:
- **Deterministic event allocation** in frequency mode: Precise control over event proportions
- **Stochastic diversity** in probability mode: Randomness while maintaining expected frequencies
- **Resource feasibility**: All generated transcripts respect exon count constraints
- **Coordinate consistency**: Genomic and transcriptomic coordinates are mutually consistent

This provides a rigorous foundation for generating gold standard datasets for AS tool evaluation, with precise control over event distributions and technical parameters.
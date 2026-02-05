# ASimulatoR Mathematical Representation

**Date:** 2026-02-05  
**Focus:** Bioinformatic science details of the ASimulatoR pipeline up to transcript construction and count matrix creation

## 1. Pipeline Overview

ASimulatoR is a splice-aware RNA-Seq data simulator that generates gold standard datasets for evaluating alternative splicing (AS) analysis tools. The pipeline consists of three main stages:

1. **Exon Superset Creation**: Merges all annotated exons from genome annotation into unspliced templates
2. **Transcript Variant Creation**: Generates transcript variants with user-defined AS event distributions
3. **RNA-Seq Read Simulation**: Simulates RNA-Seq reads using a modified version of polyester

## 2. Mathematical Notation

### 2.1 Basic Definitions

- $G = \{g_1, g_2, \ldots, g_n\}$: Set of genes in the genome annotation
- For each gene $g_i$:
  - $T_i = \{t_{i1}, t_{i2}, \ldots, t_{ik_i}\}$: Set of transcript variants
  - $E_i = \{e_{i1}, e_{i2}, \ldots, e_{im_i}\}$: Set of exons from all transcript variants (superset)
  - Each exon $e_{ij} = (chr, start_{ij}, end_{ij}, strand)$: Genomic coordinates
  - $m_i$: Number of exons in exon superset of gene $g_i$

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

**When `probs_as_freq = TRUE` and `max_genes = NULL`**, the algorithm uses all available compatible exon supersets with deterministic event assignment based on relative frequencies.

#### 4.2.1 Initial Setup

1. **Gene Length Distribution**: For all genes $G$, compute:
   $$L = \{l_i = |S_i| \mid g_i \in G\}$$
   where $l_i$ is the number of exons in gene $g_i$'s superset.

2. **Cumulative Distribution**: Create cumulative counts of genes with at least $x$ exons:
   $$C_x = \sum_{i=1}^n \mathbb{I}(l_i \geq x)$$
   This is computed as `cumsum(rev(table(factor(gene_lengths, levels = 1:max(gene_lengths)))))` - a reverse cumulative sum of the frequency table.

3. **Available Genes**: Set $N = \sum_{i=1}^n \mathbb{I}(l_i > 1)$ (genes with >1 exon, since single-exon genes cannot have AS events).

#### 4.2.2 Deterministic Event Assignment (Relative Frequency Mode)

For relative frequencies $p_c$ (where $\sum_c p_c \leq 1$):

1. **Allocate Genes per Event**: For each event combination $c$:
   $$n_c = \lfloor p_c \cdot N \rfloor$$
   Genes are deterministically assigned to events based on these counts.

2. **Create Assignment Matrix**: Build binary matrix $A$ of size $N \times |\mathcal{E}'|$ where:
   - First $n_{c_1}$ genes get event $c_1$
   - Next $n_{c_2}$ genes get event $c_2$
   - etc.
   - Remaining genes get no event (if $\sum_c p_c < 1$)

   In code: `construct_all <- sapply(names(event_probs), function(event) { rep(names(event_probs) == event, floor(event_probs[[event]] * nr_genes)) })`

3. **Flatten to List**: Convert matrix to list of event assignments per gene.

#### 4.2.3 Resource Constraint Checking

For each event combination $c$, compute minimum exon requirement:
$$r_c = \sum_{e \in c} r_e - (|c| - 1)$$

Check if sufficient genes meet this requirement:
$$\text{If } C_{r_c} < n_c \text{ for any } c$$

Where $C_{r_c}$ is the number of genes with at least $r_c$ exons.

#### 4.2.4 Iterative Adjustment

If insufficient genes:
1. Compute missing genes: $m_c = n_c - C_{r_c}$
2. Reduce total genes: $N \leftarrow N - \lceil \frac{N}{n_c} \cdot m_c \rceil$
3. Recompute allocations and check again

**Loop continues until all event requirements can be satisfied with available gene resources.**

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

The `.assign_events()` function implements a **divide-and-conquer recursive algorithm** that partitions exon indices among competing AS events:

**Input**:
- $V = \{v_1, v_2, \ldots, v_m\}$: Available exon indices (1-based positions in template)
- $R = \{(e_1, r_1), (e_2, r_2), \ldots, (e_k, r_k)\}$: Events with their minimum exon requirements $r_e$

**Algorithm**:

1. **Base Case**: If $R = \emptyset$, return $\emptyset$

2. **Recursive Step** (attempts up to `max_attempts = 10` times):
   a. **Random Event Selection**: Choose event $e \sim \text{Uniform}(R)$
   
   b. **Random Position Selection**: If attempt ≤ max_attempts:
      $$s \sim \text{Uniform}\{1, 2, \ldots, |V| - r_e + 1\}$$
      Otherwise (greedy fallback): $s = 1$
   
   c. **Assign Exons**: $f(e) = \{v_s, v_{s+1}, \ldots, v_{s+r_e-1}\}$
   
   d. **Split Remaining Exons**:
      $$V_1 = \{v_1, \ldots, v_{s-1}\}$$
      $$V_2 = \{v_{s+r_e}, \ldots, v_{|V|}\}$$
   
   e. **Partition Remaining Events**: Randomly split $R \setminus \{e\}$ between $V_1$ and $V_2$:
      $$\text{split} \sim \text{Binomial}(|R|-1, 0.5)$$
      Events assigned to $V_1$ and $V_2$ must satisfy:
      $$\sum_{e \in R_1} r_e - (|R_1| - 1) \leq |V_1|$$
      $$\sum_{e \in R_2} r_e - (|R_2| - 1) \leq |V_2|$$
   
   f. **Recursive Calls**: 
      $$f_1 = \text{assign\_events}(V_1, R_1)$$
      $$f_2 = \text{assign\_events}(V_2, R_2)$$
   
   g. **Success Condition**: If $|f_1| + |f_2| = |R| - 1$, return $\{f(e)\} \cup f_1 \cup f_2$

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
2. **Multiple Exon Skipping (mes)**: Skip intermediate exons
   $$I_{f(mes)_k} \leftarrow 0 \quad \forall k \in \{2,\ldots,|f(mes)|-1\}$$
3. **Alternative First Exon (afe)**: Swap first two exons
   $$I_1, I_2 \leftarrow I_2, I_1$$
4. **Alternative Last Exon (ale)**: Swap last two exons
   $$I_{m-1}, I_m \leftarrow I_m, I_{m-1}$$
5. **Mutually Exclusive Exons (mee)**: Choose one of two middle exons
   $$I_{f(mee)_2}, I_{f(mee)_3} \leftarrow I_{f(mee)_3}, I_{f(mee)_2}$$
6. **Intron Retention (ir)**: Merge two exons (handled separately)
7. **Alternative Splice Sites**: Modify boundaries (handled separately):
   - For $a3$: $end_{f(a3)_2} \leftarrow new\_end \in (start_{f(a3)_2} + 1, end_{f(a3)_2} - 1)$
   - For $a5$: $start_{f(a5)_1} \leftarrow new\_start \in (start_{f(a5)_1} + 1, end_{f(a5)_1} - 1)$

#### 4.4.2 Variant Transcript

The variant transcript is:
$$V_i = \{S_{ij} \mid I_j = 1\} \quad \text{with modified boundaries for A3SS/A5SS/IR}$$

### 4.5 Event Annotation Generation

The `get_event_annotation()` function creates records for each event:

1. **Exon Skipping (es)**:
   - Genomic: $(\text{start}_{skipped}, \text{end}_{skipped})$
   - Transcriptomic: $(\text{tr\_start}_{skipped}, \text{tr\_end}_{skipped})$

2. **Multiple Exon Skipping (mes)**:
   - Genomic: $(\text{start}_{skipped_1},\text{end}_{skipped_1}),\ldots,(\text{start}_{skipped_k},\text{end}_{skipped_k})$
   - Transcriptomic: $(\text{tr\_start}_{skipped_1},\text{tr\_end}_{skipped_1}),\ldots$

3. **Intron Retention (ir)**:
   - Genomic: $(\text{end}_{upstream} + 1, \text{start}_{downstream} - 1)$
   - Transcriptomic: $(tr\_start_{retained}, tr\_end_{retained})$ where:
     $$tr\_start_{retained} = tr\_start_{upstream} + \Delta$$
     $$\Delta = \begin{cases}
     \text{end}_{upstream} - \text{end}_{ri} & \text{if strand} = - \\
     \text{start}_{ri} - \text{start}_{upstream} & \text{otherwise}
     \end{cases}$$

4. **Alternative Splice Sites**:
   - For $a3$: $C_g = (\min(new\_end, old\_end) + 1, \max(new\_end, old\_end))$
     $$C_t = (tr\_start_{exon}, tr\_start_{exon} + |new\_end - old\_end| - 1)$$
   - For $a5$: $C_g = (\min(new\_start, old\_start), \max(new\_start, old\_start) - 1)$
     $$C_t = (tr\_end_{exon} - |new\_start - old\_start| + 1, tr\_end_{exon})$$

5. **Mutually Exclusive Exons (mee)**: Annotates which exon is included in variant vs template

6. **Alternative First/Last Exons (afe/ale)**: Annotates swapped exons

## 5. RNA-Seq Read Simulation (Overview)

### 5.1 Polyester Integration

ASimulatoR uses a modified version of polyester with parameters:

- **Sequencing depth**: $D$ reads per sample
- **Read length**: $L$ bases
- **Error rate**: $\epsilon$
- **Fragment length distribution**: $FL \sim \text{Normal}(\mu_f, \sigma_f^2)$
- **Strand-specific**: Boolean

### 5.2 Read Generation Model

For each transcript $t$ with length $l_t$:

1. **Expected reads**: $r_t = D \cdot \frac{l_t}{\sum_{t'} l_{t'}} \cdot \text{expression\_factor}$

2. **Fragment generation**: Sample $r_t$ fragment start positions uniformly from $[1, l_t - FL + 1]$

3. **Read simulation**: For each fragment, generate paired-end reads with:
   - Sequencing errors: Each base error with probability $\epsilon$
   - Quality scores: Phred-scaled based on error probability

### 5.3 Technical Biases

1. **PCR Duplicates**: Fraction $d$ of reads duplicated, duplication count $\sim \text{Poisson}(\lambda=1)$

2. **Adapter Contamination**: For fragments shorter than $L$, append adapter sequence

## 6. Mathematical Summary

### 6.1 Key Equations

1. **Exon superset size**: $|S_i| = |\text{reduce}(\bigcup_{t \in T_i} E_{it})|$

2. **Minimum exons for event combination**: 
   $$\min\_exons(E) = \sum_{e \in E} r_e - (|E| - 1)$$
   where $r_e$ is minimum exon requirement for event $e$

3. **Event assignment probability**: For probability mode:
   $$P(\text{gene } g_i \text{ has event } e) = p_e \cdot \mathbb{I}(|S_i| \geq r_e)$$

4. **Variant construction**: 
   $$V_i = \text{apply\_events}(S_i, f, I)$$
   where $f$ is event-exon mapping and $I$ is inclusion vector

### 6.2 Algorithm Complexity

- **Exon superset creation**: $O(\sum_i |T_i| \cdot \log |E_i|)$
- **Event assignment**: $O(N \cdot |\mathcal{E}'|)$
- **Variant construction**: $O(\sum_i |S_i| \cdot |E_i|)$

## 7. Biological Interpretation

The mathematical model captures:

1. **Exon combinatorics**: All possible exon combinations from annotated transcripts
2. **AS event constraints**: Biological feasibility of event types based on exon count
3. **Strand awareness**: Proper orientation of negative strand genes
4. **Coordinate systems**: Both genomic and transcriptomic coordinates for AS analysis

This provides a rigorous foundation for generating gold standard datasets for AS tool evaluation, with precise control over event distributions and technical parameters.
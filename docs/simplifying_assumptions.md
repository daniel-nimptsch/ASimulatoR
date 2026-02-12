1. Random uniform selection for gene assignment - Genes are selected uniformly at random from the pool of compatible genes
2. Random event placement - Events are placed randomly on exons within each gene
3. Fixed fold changes for differential expression - Only uses 2x fold change, fixed at 2
4. Uniform splice site position for A3SS/A5SS - Alternative splice sites are uniformly distributed within exon boundaries
5. Uniform exon selection for afe/ale/mee - Exons are selected uniformly at random
6. Fixed 50% probability for isoform switches - When multiple transcripts exist
7. Length-based expression model only - No biological variation in expression beyond transcript length
8. No GC bias simulation - Explicitly disabled (frag_GC_bias = 'none')
9. No sequence-specific effects - Splice site strength, exon definition, etc. not modeled
10. Contiguous event blocks only - Events must be on contiguous exons
11. Non-overlapping events - Events cannot share exons (unless multi_events_per_exon enabled)
12. Simplified intron retention - Just merges two exons, doesn't model complex retention patterns
13. Fixed PCR duplication distribution - Poisson(λ=1) for number of duplications
14. Uniform fragment distribution - Fragments uniformly distributed across transcripts
15. No read quality variation beyond error rate - All reads have same expected quality based on error rate
16. Fixed maximum attempts for event placement - 10 attempts before greedy fallback
17. Simplified biological realism - Events assigned randomly without considering splice site motifs, RBP binding sites, etc.
18. No tissue-specific or condition-specific expression - Expression determined purely by length and fold changes
19. No evolutionary constraints - Doesn't consider whether simulated events are biologically plausible
20. Fixed splice site modification range - A3SS/A5SS positions sampled uniformly from exon interior (excluding endpoints)

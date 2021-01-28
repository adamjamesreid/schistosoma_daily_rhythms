# Files relating to Rawlinson et al. Daily rhythms in the adult transcriptomes of the blood fluke, Schistosoma mansoni

### schisto_cluster_enrichment.py
This Python script was used to determine enrichment of genes with cycling transcripts amongst marker genes from single-cell RNA-seq data cell types (Wendt et al., 2020; https://science.sciencemag.org/content/369/6511/1644.abstract)

It is run in this way:

```
python schisto_cluster_enrichment.py <markers> <cycling list>
```
e.g.
```
python schisto_cluster_enrichment.py rmv27_50_markers_RN.csv pooled_male_cycling_genes.txt
```

### rmv27_50_markers_RN.csv
List of marker genes for different cell clusters from Wendt et al. (2020)

### Files of cycling genes
pooled_male_cycling_genes.txt
pooled_female_cycling_genes.txt
pooled_head_cycling_genes.txt

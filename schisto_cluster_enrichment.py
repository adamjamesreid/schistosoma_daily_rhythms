# Determine which schisto single-cell clusters are enriched for cycling genes
# Input: gene:cluster, cycling genes

import sys
import scipy, scipy.stats
from scipy.stats import hypergeom
import statsmodels.stats.multitest

fdr_cutoff = 0.05
logfc_threshold = 0.25 # marker fc threshold for inclusion
min_pct = 0.25 # marker pct threshold for inclusion

clusters = sys.argv[1]
cycling = sys.argv[2]

# get gene:cluster in dict
cl = open(clusters)
clusters = dict()
gene_count = 0
for x in cl.readlines():
    x = x.rstrip()
    v = x.split(',')
    if v[0] == 'gene':
        continue
    v[6] = v[6].replace('"', '')
    v[0] = v[0].replace('"', '')
    if len(v) < 8:
        continue
    # Exclude markers not passing thresholds
    if float(v[3]) < logfc_threshold:
        continue
    if float(v[4]) < min_pct:
        continue
    if v[7] not in clusters:
        clusters[v[7]] = list()
    clusters[v[7]].append(v[0])
    gene_count = gene_count + 1
    #print(v[6], v[0])

# get cycling genes in list
cy = open(cycling)
cycling = list()
for x in cy.readlines():
    x = x.rstrip()
    v = x.split('\t')
    if v[0] != 'id':
        v[0] = v[0].replace('_', '-')
        cycling.append(v[0])

# Loop through clusters and get num genes in cluster, num cycling in cluster, total genes, total cycling
pvals = list()
results = list()
cycling_markers = dict()
for c in clusters:
    #print(c)
    genes_in_cluster = len(clusters[c]) # N
    cycling_in_cluster = len(set(clusters[c]).intersection(set(cycling)))
    cycling_markers[c] = list()
    cycling_markers[c] = set(clusters[c]).intersection(set(cycling))
    total_cycling = len(cycling) # n
    if cycling_in_cluster < 1:
        continue
    dens = hypergeom.cdf(cycling_in_cluster, gene_count, total_cycling, genes_in_cluster)
    p = 1 - dens
    pvals.append(p)
    #print(c, cycling_in_cluster, gene_count, total_cycling, genes_in_cluster, p)
    results.append([c, cycling_in_cluster, gene_count, total_cycling, genes_in_cluster, p])

if len(pvals) > 0:
  qvals = statsmodels.stats.multitest.multipletests(pvals, alpha = 0.05, method='fdr_bh')
  for x in range(0, len(results)):
    cluster_str = results[x][0]
    results[x][0] = '"' + cluster_str + '"'
    res_str = '\t'.join([str(s) for s in results[x]])
    genes = ' '.join(cycling_markers[cluster_str])
    fdr = qvals[1][x]
    if (fdr <= fdr_cutoff):
      print(res_str, fdr, genes, sep='\t')

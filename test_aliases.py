"""Quick check: what aliases do Cluster 0 markers have?"""
import sys
sys.path.insert(0, '/home/kwy7605/data_61/scRNAseq-tools/src/scrnaseq_tools/templates/claude/skills/Annotation-agent/tools')
from dynamic_knowledge import _get_cached_aliases

for gene in ['FCER2', 'BTLA', 'BANK1', 'CD200', 'LTB', 'KLF2', 'FOXO1', 'BACH2']:
    aliases = _get_cached_aliases(gene)
    print(f"{gene}: {aliases}")

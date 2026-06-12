# Processed Data Tables

These CSV files are processed, non-sample-level source tables used by the plotting scripts.

## Files

- `i4_fgsea_hallmark.csv`: Inspiration4 fGSEA results for MSigDB Hallmark gene sets.
- `i4_fgsea_c4.csv`: Inspiration4 fGSEA results for MSigDB C4 computational cancer gene sets.
- `i4_fgsea_c6.csv`: Inspiration4 fGSEA results for MSigDB C6 oncogenic signature gene sets.
- `i4_key_cancer_genes_log2fc.csv`: Key cancer gene log2FC and p-value summary across PBMC and immune cell subsets.

## Common fGSEA Columns

- `pathway`: MSigDB pathway or gene-set name.
- `pval`: Nominal enrichment p-value.
- `padj`: Multiple-testing adjusted p-value.
- `ES`: Enrichment score.
- `NES`: Normalized enrichment score.
- `size`: Gene-set size used by fGSEA.
- `leadingEdge`: Semicolon-delimited leading-edge genes.
- `timepoint`: Post-flight comparison.
- `celltype`: PBMC or immune cell subset.
- `info`: Gene-set collection.

## Key Cancer Gene Columns

- `Gene`: Gene symbol.
- `Cell_Type`: PBMC or immune cell subset.
- `log2FC_R1_vs_Preflight`: Immediate post-flight log2 fold change versus pre-flight.
- `log2FC_R45_R82_vs_Preflight`: Long-term post-flight log2 fold change versus pre-flight.
- `pval_*` and `padj_*`: Nominal and adjusted p-values for each comparison.
- `In_OncoKB`, `In_HRR`, `In_NHEJ`, `In_Key_Cancer`: Gene-set membership flags used for manuscript panels.

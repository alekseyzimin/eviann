# EviAnn
An evidence-based eukaryotic annotation pipeline.  This is the eviann submodule.  For the EviAnn releases please go to EviAnn_release repository.

## External CDSs

By default, external CDSs supplied with `-c` (or `--cds`) are trusted and bypass the splice score thresholds when they do not match assembled transcripts.

Use `--untrusted-cds` with `-c` to apply the same Markov/WAM thresholds used for unmatched protein-derived candidates. Both splice scores must exceed the existing threshold. This only affects external CDSs unmatched to transcripts; downstream external CDS handling is unchanged. Single-exon models have no splice junctions and pass this filter.

Example:

```bash
eviann.sh -g genome.fa -e ests.fa -c external_cds.gff --untrusted-cds
```

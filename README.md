# EviAnn
An evidence-based eukaryotic annotation pipeline.  This is the eviann submodule.  For the EviAnn releases please go to EviAnn_release repository.

## External CDSs

By default, external CDSs supplied with `-c` (or `--cds`) are treated as trusted evidence.

For external CDSs that are not high-confidence, add `--untrusted-cds` to apply a subset of EviAnn's splice-site filters. This filtering applies only to CDSs without matching transcript evidence. Single-exon models have no splice sites and pass this filter.

Example:

```bash
eviann.sh -g genome.fa -e ests.fa -c external_cds.gff --untrusted-cds
```

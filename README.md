# Hopper et al. 2016 — Grapevine Leaf Dehydration Transcriptomics

## Publication

**Title:** Transcriptomic network analyses of leaf dehydration responses identify highly connected ABA and ethylene signaling hubs in three grapevine species differing in drought tolerance

**Authors:** Daniel W. Hopper, Ryan Ghan, Karen A. Schlauch, Grant R. Cramer

**Journal:** *BMC Plant Biology*, 2016, 16(1):118

**DOI:** [10.1186/s12870-016-0804-6](https://doi.org/10.1186/s12870-016-0804-6)

**PMID:** [27215785](https://pubmed.ncbi.nlm.nih.gov/27215785/)

---

## Study Overview

This study used time-series transcriptomics and weighted gene co-expression network analysis (WGCNA) to characterize how grapevine leaves respond to dehydration stress across three *Vitis* species that differ in drought tolerance:

| Genotype | Drought Tolerance |
|---|---|
| *Vitis champinii* cv. Ramsey | High (drought-tolerant) |
| *Vitis vinifera* cv. Cabernet Sauvignon | Intermediate |
| *Vitis riparia* cv. Riparia Gloire | Low (drought-sensitive) |

Leaves were sampled at multiple time points across a controlled 24-hour dehydration assay. Transcriptomic profiling was performed using a *Vitis vinifera* microarray platform.

---

## Key Findings

- **~11,000 transcripts** changed significantly with the genotype × treatment interaction (non-parametric 3-way ANOVA)
- **~18,237 transcripts** responded significantly to treatment; **~23,656** to time
- **WGCNA** identified **30 co-expression modules** with distinct biological enrichments
- Key enriched pathways: photosynthesis, phenylpropanoid metabolism, ABA signaling, ethylene signaling
- **VviABI5** and **VviABF2** were identified as highly connected hubs in ABA and ethylene signaling modules
- **VviABI5** showed an early, strong response in drought-tolerant Ramsey with minimal response in drought-sensitive Riparia Gloire
- **VviSnRK1** and related sugar/hormone signaling genes were co-hub candidates
- Results support substantial crosstalk between ABA and ethylene pathways during dehydration stress

---

## Repository Contents

### `Hopper16.pdf`
The full published article (Hopper et al., 2016, *BMC Plant Biology*).

### `Hopper16Supp3.xlsx`
Supplementary data file 3. Contains the combined leaf gene expression dataset analyzed with non-parametric 3-way statistics.

**Size:** ~20.4 MB

#### Sheet 1: `Hopper.Combined.Leaf.NonParam3w`

Combined leaf expression dataset across all three grapevine genotypes, analyzed with a non-parametric 3-way statistical approach.

| Column(s) | Content |
|---|---|
| A | Probe/transcript ID |
| B–E | Additional identifiers / integer indices |
| F–CY (98 columns) | Normalized expression values (log-scale, ~6–14 range) |

- **Rows:** 29,551 (one per probe/transcript)
- **Samples:** 98 expression columns spanning all genotypes, treatment conditions, and time points

#### Sheet 2: `Sheet1`

Summary/reference sheet with a small subset (4 rows × 47 columns) containing:
- Sample reference IDs
- Expression measurements for representative transcripts
- Statistical probability values (p-values and model fit scores) in columns AH–AU

---

## Methods Summary

- **Experimental design:** Leaf detachment dehydration assay; three grapevine species; multiple time points over 24 hours
- **Platform:** *Vitis vinifera* microarray
- **Statistical analysis:** Non-parametric 3-way ANOVA (genotype × treatment × time)
- **Network analysis:** WGCNA (Weighted Gene Co-expression Network Analysis)
- **GO enrichment:** Overrepresentation analysis of gene ontology categories per WGCNA module

---

## Keywords

ABA · ABI5 · Dehydration · Drought tolerance · Ethylene · Grapevine · Network analysis · Transcriptomics · *Vitis* · WGCNA

---

## Citation

```
Hopper DW, Ghan R, Schlauch KA, Cramer GR.
Transcriptomic network analyses of leaf dehydration responses identify highly
connected ABA and ethylene signaling hubs in three grapevine species differing
in drought tolerance.
BMC Plant Biology. 2016;16(1):118.
doi:10.1186/s12870-016-0804-6
```

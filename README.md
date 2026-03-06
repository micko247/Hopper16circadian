# Hopper et al. 2016 — Grapevine Leaf Dehydration Transcriptomics

## Publication

**Title:** Transcriptomic network analyses of leaf dehydration responses identify highly connected ABA and ethylene signaling hubs in three grapevine species differing in drought tolerance

**Authors:** Daniel W. Hopper, Ryan Ghan, Karen A. Schlauch, Grant R. Cramer

**Affiliation:** Department of Biochemistry and Molecular Biology, University of Nevada, Reno, NV 89557, USA

**Journal:** *BMC Plant Biology*, 2016, 16(1):118

**DOI:** [10.1186/s12870-016-0804-6](https://doi.org/10.1186/s12870-016-0804-6) | **PMID:** [27215785](https://pubmed.ncbi.nlm.nih.gov/27215785/) | **GEO:** [GSE78920](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE78920)

---

## Study Overview

This study used time-series transcriptomics and weighted gene co-expression network analysis (WGCNA) to characterize how grapevine leaves respond to dehydration stress across three *Vitis* species differing in drought tolerance.

### Grapevine Genotypes

| Genotype | Species | Drought Tolerance |
|---|---|---|
| Ramsey | *Vitis champinii* Planch. | High — native to hot, dry regions of Texas |
| Cabernet Sauvignon | *Vitis vinifera* L. cv. Cabernet Sauvignon clone 8 | Intermediate |
| Riparia Gloire | *Vitis riparia* Michx. | Low — native to wet riparian areas |

---

## Experimental Design

- **Design:** 2 × 2 × 5 factorial — Genotype (3) × Treatment (control / dehydration) × Time (1, 2, 4, 8, 24 h)
- **Assay:** Fully developed leaves excised and placed in a sealed dehydration box with NaCl solution to control humidity; frozen in liquid nitrogen at each time point
- **Controls:** Control leaves were harvested from the same plant at the corresponding time of day as the stress leaves, explicitly to account for **circadian effects on transcript abundance**
- **Replication:** n = 3 biological replicates per treatment × time point combination (90 arrays total; 3 excluded as QC outliers → 87 final arrays)
- **Platform:** NimbleGen Grape Whole-Genome Microarray (090818 Vitis exp HX12)
- **RNA extraction:** CTAB-based method; quality assessed by Nanodrop and Agilent Bioanalyzer
- **Greenhouse conditions:** 16 h/8 h light/dark, 28 °C/18 °C, ≥400 μE m⁻² s⁻¹

---

## Statistical Methods

| Method | Details |
|---|---|
| ANOVA | Non-parametric extension of Kruskal-Wallis rank sum test on log-transformed, normalized data |
| FDR correction | Benjamini-Hochberg; adjusted p ≤ 0.05 |
| PCA | Covariance matrix |
| WGCNA | `power=16`, `type="signed hybrid"`, `corFnc="bicor"`, `minClusterSize=30` |
| GO enrichment | BiNGO (Cytoscape); hypergeometric test; Benjamini-Hochberg FDR p ≤ 0.05 |
| Ortholog mapping | Gramene release 44 (January 2015) |

---

## Key Findings

### Scale of Transcriptional Response (adjusted p ≤ 0.05)

| ANOVA Effect | Significant Transcripts |
|---|---|
| Genotype (G) | 28,030 |
| Time | 23,656 |
| Genotype × Time | 24,543 |
| Treatment (TRT) | 18,237 |
| Treatment × Time | 17,488 |
| **Genotype × Treatment** | **11,436** |
| Genotype × Treatment × Time | 6,285 |

### ABA Signalling

- **VviNCED6** (ABA biosynthesis): increased within 1 h in North American species (Ramsey, Riparia Gloire) but not Cabernet Sauvignon
- **VviCYP707A4** (ABA catabolism): decreased log₂ fold > 4 at 4 h in all genotypes
- **VviABCG25 / VviABCG40** (ABA transporters): both increased; VviABCG40 increased ~5-fold in Riparia Gloire within 1 h
- **VviABI5**: most discriminating gene — dramatic early increase in Ramsey, minimal response in Riparia Gloire
- **VviABF2, VviOST1, VviHAI1, VviAHG3**: all differentially expressed between genotypes

### Ethylene Signalling

- Multiple **ACS** genes induced (VviACS2, 4, 6, 7, 8-like); VviACS2 highest in Riparia Gloire at 1 h
- **WRKY transcription factors** VviWRKY33 and VviWRKY40 increased rapidly within 1 h in all genotypes
- 91 of 130 AP2/ERF superfamily members changed significantly across genotypes
- Riparia Gloire showed a rapid, heightened ethylene response; Ramsey showed a suppressed ethylene response

### WGCNA Network Modules

WGCNA identified **30 co-expression modules**. Two modules with highest relevance to drought tolerance:

**Yellow3 module** (enriched: gas/oxygen transport, p = 1.28 × 10⁻⁴)
- Top hub: VIT_10s0071g00840 (TPR domain protein)
- Key hubs (kME > 0.80): VviABI5, VviABF2, VviSnRK1, VviCAD, VviBCH, VviGDSL1, VviUSP
- 81 hub genes with kME > 0.80

**Lightsteelblue module** (enriched: ethylene signalling GO)
- Key hubs: VviABI5, VviABF2, VviSnRK1, VviERF1, VviNCED3, VviRAV2, VviGID1B, VviRAN1, VviSTP3, VviIDD2

VviABI5 and VviABF2 are highly connected hubs in **both** modules, providing network-level evidence for ABA–ethylene crosstalk during dehydration stress.

---

## Repository Contents

### `Hopper16.pdf`
Full published article (Hopper et al., 2016, *BMC Plant Biology*).

### `Hopper16Supp3.xlsx`
Supplementary data file 3 (~20.4 MB). Contains the combined leaf gene expression dataset analyzed with non-parametric 3-way statistics.

#### Sheet 1: `Hopper.Combined.Leaf.NonParam3w`

Combined leaf expression data across all three genotypes, analyzed with non-parametric 3-way ANOVA (genotype × treatment × time).

| Column(s) | Content |
|---|---|
| A | Probe / transcript ID |
| B–E | Additional integer identifiers / indices |
| F–CY (98 columns) | Normalized, log-transformed expression values (~6–14 range) |

- **Rows:** 29,551 (one per probe/transcript)
- **Samples:** 98 expression columns covering all genotypes, treatment conditions, and time points

#### Sheet 2: `Sheet1`

Small summary/reference sheet (4 rows × 47 columns) containing sample reference IDs, representative expression measurements, and statistical probability values (p-values and model fit scores in columns AH–AU).

---

## Keywords

ABA · ABI5 · Dehydration · Drought tolerance · Ethylene · Grapevine · Network analysis · Transcriptomics · *Vitis* · WGCNA

---

## Funding

- USDA-NIFA Hatch Grant NEV00345
- NIH NIGMS P20GM103440

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

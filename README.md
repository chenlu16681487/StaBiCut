# StaBiCut
**StaBiCut** (Stable and Biologically aware Cutoff prioritization) is an R-based framework for **stability- and biology-aware survival cutoff prioritization** in candidate biomarker studies.

<img width="523" height="208.3" alt="Stabicut_logo" src="https://github.com/user-attachments/assets/cff56554-0402-4642-bf6e-a7197daf5df1" />




Rather than proposing a new survival test, StaBiCut is designed as a **decision layer on top of classical survival analysis**. It helps users determine whether a survival-associated expression cutoff is not only statistically separative, but also **directionally coherent, distributionally plausible, locally reproducible under resampling, and practically interpretable**.

This repository contains the code used for the current StaBiCut implementation and the analyses underlying its application to the TCGA colorectal cancer cohort (COAD and READ).

---
  
  ## Overview
  
  Classical survival workflows can identify a cutoff that yields an apparently strong survival split in one cohort. However, such a cutoff may still be:
  
- statistically fragile under resampling,
- driven by expression tails or highly imbalanced partitions,
- inconsistent with prior biological or tumor–normal directional evidence, or
- difficult to prioritize across multiple candidate genes.

StaBiCut was developed to address this practical gap. It retains the standard survival-analysis backbone (`coxph`, `survdiff`, Kaplan–Meier estimation), while adding a structured prioritization framework based on:
  
1. **direction-aware candidate filtering**,  
2. **deterministic cutpoint scanning**,  
3. **local bootstrap re-support of the selected cutoff**,  
4. **post-selection plausibility and consistency metrics**, and  
5. **fixed-weight composite ranking across genes**.

Accordingly, StaBiCut should be interpreted as a **prior-aware, direction-constrained cutoff prioritization framework**, rather than as a purely data-driven search for the most significant threshold.

---
  
  ## Conceptual positioning
  
  StaBiCut does **not** replace `survival`, `survminer`, or `maxstat`. Instead, it operates one layer above them.

- Classical tools answer: **Can a model be fitted? Can a cutoff be found?**
  - StaBiCut answers: **Is this cutoff worth trusting and prioritizing?**
  
  In this sense, StaBiCut is best understood as a framework for **stability-aware biomarker cutoff prioritization** in translational studies.

---
  
  ## Core design principles
  
  StaBiCut was built around the following principles:
  
- **Prior-aware**: when available, expected hazard direction can be constrained using gene-specific prior information.
- **Deterministic**: cutoff selection follows a fixed hierarchical rule, avoiding arbitrary tie resolution.
- **Stability-audited**: bootstrap is used to test whether the originally selected cutoff and its local neighborhood are reproducibly supported.
- **Not tail-driven**: cutoff plausibility and group balance are explicitly evaluated after selection.
- **Conflict-aware**: tumor–normal trend and survival direction are jointly considered.
- **Ranking-oriented**: outputs are designed to support prioritization across multiple candidate genes, not merely single-gene significance claims.

---
  
  ## Workflow overview
  
  A workflow overview is provided here:
  
 <img width="2851" height="1613" alt="StaBiCut_workflow" src="https://github.com/user-attachments/assets/c4ab7778-48ab-4162-8ed5-f65fb04cf576" />


For the mathematical definitions of the scoring components and stability metrics, please see:
  
  - `Supplementary Methods`
- future vignette / method note linked in this repository

---
  
  ## Quick start in 3 commands
  
  ```r
source("R/modules_StaBiCut_v2.R")
source("R/run_StaBiCut_v2.R")
source("R/Panel_helper_StaBiCut_v2.R")

Then run:
set.seed(1)
res <- run_batch_sur_cutpoint_analysis_v2(
  exprset = mrna_expr_tpm,
  geneset = geneset,
  clin = clinicalSE,
  gene_prior_table = gene_prior_table,
  force_direction = FALSE,
  n_boot = 1000,
  minprop = 0.25,
  save_boot_rds = TRUE,
  boot_dir = "./output/boot_rds",
  save_plots = TRUE,
  plot_dir = "./output",
  seed = 1)  

And for cross-seed robustness:
all_runs <- run_multiseed_stability_v2(
    seeds = 1:20,
    n_boot = 1000,
    exprset = mrna_expr_tpm,
    geneset = geneset,
    clin = clinicalSE,
    gene_prior_table = gene_prior_table,
    force_direction = FALSE,
    minprop = 0.25,
    save_each_seed_rds = TRUE,
    out_dir = "./output/stability_n1000_with_rds") 
```

Method summary
In the current implementation, StaBiCut operates as follows:
1.	Expression preprocessing
Expression matrices are transformed to log2(TPM + 1). Genes with extremely low mean expression are removed. Survival analyses are restricted to tumor samples. Tumor expression values are winsorized at the 5th and 95th percentiles before cutoff scanning. 
2.	Direction-prior layer
If a valid gene-specific prior table is provided, that expected direction is used directly. Otherwise, direction can be inferred only when tumor–normal expression shift and the sign of the continuous Cox coefficient are concordant. If reliable directional evidence is unavailable, no direction constraint is imposed. 
3.	Deterministic cutoff scan
Candidate cutoffs are scanned within the central expression range defined by minprop (default: 0.25). Candidates that yield highly imbalanced groups, insufficient events, failed Cox fitting, or unstable regression estimates are discarded. When direction is specified, only direction-compatible cutoffs are retained. 
4.	Fixed cutoff estimation
The final cutoff is selected by a deterministic hierarchical rule that prioritizes the Cox Wald statistic, then the corresponding P value, and then the distance of HR from 1. The retained cutoff is used for dichotomized Cox modeling and log-rank testing. 
5.	Local bootstrap re-support
Bootstrap resampling is used to evaluate whether the selected cutoff itself is reproducibly supported under perturbation. Rather than re-optimizing a new global optimum in each bootstrap replicate, StaBiCut rescans candidate cutoffs only within a local window centered on the originally selected cutoff. 
6.	Five-component scoring and prioritization
Each gene is summarized by five components:
o	bootstrap stability, 
o	cutoff plausibility, 
o	bootstrap hazard-direction consistency, 
o	tumor–normal directional concordance, 
o	group balance. 
These are integrated into a weighted composite score for cross-gene prioritization.
7.	Cross-seed robustness analysis
The full pipeline can be repeated across multiple random seeds. Cross-run stability is then summarized by dispersion of composite scores, dispersion of selected cutoffs, and TopK ranking frequency. 

Scoring components
In the current implementation, the composite prioritization score combines:
- Bootstrap stability (weight 0.40) 
- Cutoff plausibility (weight 0.20) 
- Bootstrap hazard-direction consistency (weight 0.15) 
- Tumor–normal directional concordance (weight 0.15) 
- Group balance (weight 0.10) 
Higher scores indicate cutoffs that are more reproducible, more biologically coherent, more distributionally plausible, and more robust and practically informative for downstream interpretation.
Mathematical definitions are documented in the Supplementary Methods and may later be expanded into a repository vignette.

Repository structure
The current public release is organized as a script-style research repository, with future migration to a more package-like structure if needed.
```
StaBiCut/
├── R/
│   ├── modules_StaBiCut_v2.R
│   ├── run_StaBiCut_v2.R
│   └── Panel_helper_StaBiCut_v2.R
├── examples/
│   ├── example_run.R
│   ├── example_gene_prior_table.csv
│   └── example_output_dictionary.csv
├── figures/
│   └── StaBiCut_workflow_overview.png
├── README.md
├── LICENSE
├── CITATION.cff
├── sessionInfo.txt
└── .gitignore
```

At present, the core implementation is centered on:
- run_batch_sur_cutpoint_analysis_v2() 
- scan_cutpoints_v2() 
- bootstrap_cutoffs_v2() 
- run_multiseed_stability_v2() 
- plot_gene_panel_main_v2() 
- plot_multi_gene_sheet_main_v2() 
- plot_multiseed_composite_summary_v2() 
- export_stability_excel_v2() 

Installation
Option 1: source the scripts directly
source("R/modules_StaBiCut_v2.R")
source("R/run_StaBiCut_v2.R")
source("R/Panel_helper_StaBiCut_v2.R")

Option 2: install from GitHub
(Recommended only after the repository is formally reorganized as a standard R package.)
# remotes::install_github("chenlu16681487/StaBiCut")

Dependencies
StaBiCut currently relies on the following R packages:
- survival 
- survminer 
- ggplot2 
- patchwork 
- dplyr 
- tidyr 
- openxlsx 
- splines 
- forcats 

Exact dependency versions are not hard-coded in this README. Please record them in sessionInfo.txt when generating publication-linked results.

Input requirements
The main runner expects:
- Expression matrix on TPM scale, with genes in rows and samples in columns. 
- Clinical data frame matched to the tumor samples, containing survival information. 
- Gene set specifying the candidate genes to be evaluated. 
- Optional gene prior table with columns: 
o	Gene: Symbol 
o	expected_dir：prior expected direction of genes 
Supported direction labels are:
- protective_high 
- adverse_high 
The current implementation assumes TCGA-style sample barcodes when identifying tumor versus normal samples. For non-TCGA datasets, users should provide explicit sample-type labels or preprocess the dataset before running the current scripts.
Example object structures used in the present manuscript-linked workflow:
```
- Expression matrix on TPM scale:
> mrna_expr_tpm[1:5,1:3]
        TCGA-AA-3688-01A-01R-0905-07 TCGA-G4-6298-01A-11R-1723-07 TCGA-AA-3672-01A-01R-0905-07
MT-CO2                      25914.66                     40456.34                     15299.70
MT-CO3                      26244.67                     37331.63                     12108.39
MT-ND4                      13102.97                     25710.66                     11878.03
MT-CO1                      11398.30                     14542.93                     14488.23
MT-ATP6                     11133.31                     13070.55                     14943.64

- Clinical data frame:
> clinicalSE[1:5,1:3]
                                                  barcode      patient           sample
TCGA-AA-3688-01A-01R-0905-07 TCGA-AA-3688-01A-01R-0905-07 TCGA-AA-3688 TCGA-AA-3688-01A
TCGA-G4-6298-01A-11R-1723-07 TCGA-G4-6298-01A-11R-1723-07 TCGA-G4-6298 TCGA-G4-6298-01A
TCGA-AA-3672-01A-01R-0905-07 TCGA-AA-3672-01A-01R-0905-07 TCGA-AA-3672 TCGA-AA-3672-01A
TCGA-G4-6314-01A-11R-1723-07 TCGA-G4-6314-01A-11R-1723-07 TCGA-G4-6314 TCGA-G4-6314-01A
TCGA-A6-2682-01A-01R-1410-07 TCGA-A6-2682-01A-01R-1410-07 TCGA-A6-2682 TCGA-A6-2682-01A

- Gene set:
> geneset 
 [1] "SPOCK2"  "PYCR1"   "CA4"     "CES1"    "ABCB1"   "ZG16"    "TNXB"    "HMCN2"   "MEP1A"  
[10] "SLC37A2" "CHGB"

- Optional gene prior table:
> gene_prior_table[1:5,1:2]
    Gene    expected_dir
1 SPOCK2    adverse_high
2  PYCR1    adverse_high
3    CA4 protective_high
4   CES1 protective_high
5  ABCB1 protective_high
```

Benchmark gene set used in the current manuscript-linked workflow

The current TCGA-CRC benchmark workflow evaluates the following 11 candidate genes:

SPOCK2 / PYCR1 / CA4 / CES1 / ABCB1 / ZG16 / TNXB / HMCN2 / MEP1A / SLC37A2 / CHGB

These genes were pre-nominated from longitudinal CAC transcriptomic and proteomic analyses and were therefore analyzed with explicit gene-specific direction priors in the present study.

Minimal example
```
source("R/modules_StaBiCut_v2.R")
source("R/run_StaBiCut_v2.R")
source("R/Panel_helper_StaBiCut_v2.R")

geneset <- c("SPOCK2","PYCR1","CA4","CES1","ABCB1","ZG16",
             "TNXB","HMCN2","MEP1A","SLC37A2","CHGB")

gene_prior_table <- data.frame(
  Gene = geneset,
  expected_dir = c(rep("adverse_high", 2), rep("protective_high", 9)),
  stringsAsFactors = FALSE
)

set.seed(1)

res <- run_batch_sur_cutpoint_analysis_v2(
  exprset = mrna_expr_tpm,
  geneset = geneset,
  clin = clinicalSE,
  gene_prior_table = gene_prior_table,
  force_direction = FALSE,
  n_boot = 1000,
  minprop = 0.25,
  save_boot_rds = TRUE,
  boot_dir = "./output/boot_rds",
  save_plots = TRUE,
  plot_dir = "./output",
  seed = 1
)

head(res$results_df)

Cross-seed stability example
all_runs <- run_multiseed_stability_v2(
  seeds = 1:20,
  n_boot = 1000,
  exprset = mrna_expr_tpm,
  geneset = geneset,
  clin = clinicalSE,
  gene_prior_table = gene_prior_table,
  force_direction = FALSE,
  minprop = 0.25,
  save_each_seed_rds = TRUE,
  out_dir = "./output/stability_n1000_with_rds")
```

Main outputs
The main pipeline returns a list containing:
- results_df
Per-gene summary table including HR, confidence interval, Cox P value, cutoff, log-rank P value, bootstrap metrics, plausibility metrics, and composite score. 
- df_cache
Per-gene analysis tables used for fitted survival and plotting steps. 
- boot_cache
Per-gene bootstrap results, including retained bootstrap cutoffs and summary statistics. 
- scan_cache
Per-gene cutoff scan tables containing candidate cutoffs and scan statistics. 
- exprset_full
Full log2-transformed expression matrix. 
- exprset_tumor
Tumor-only log2-transformed expression matrix. 
The repository also supports:
- per-gene four-panel summaries, 
- multi-gene panel sheets, 
- cross-seed stability heatmaps, 
- composite score summary plots, 
- Excel export of stability statistics and column dictionaries. 
A repository-level output dictionary should be provided as either:
- a generated CSV exported from build_column_dictionary_resultsdf_v2(), or 
- a short vignette / method note describing the result columns. 

Interpretation of selected cutoffs
StaBiCut-selected cutoffs should be interpreted as practical working thresholds for cohort-level stratification under the stated decision rules, not as exact biological boundaries.
A high-ranking gene under StaBiCut is not simply one with a small P value. Rather, it is a gene whose selected cutoff is jointly supported by:
- reproducibility under local resampling, 
- acceptable placement within the tumor-expression distribution, 
- stable hazard direction, 
- coherence with tumor–normal trend, and 
- workable group balance. 

Reproducibility and benchmark configuration
The publication-linked TCGA-CRC application used:
- TCGA COAD + READ 
- 583 valid tumor samples 
- tumor-only survival analysis 
- log2(TPM + 1) expression 
- winsorization at the 5th and 95th percentiles 
- minprop = 0.25 
- n_boot = 1000 
- 20 random seeds for cross-run robustness assessment 
- gene-specific direction priors derived from the following 11 CAC multi-omics candidates:
  
SPOCK2, PYCR1, CA4, CES1, ABCB1, ZG16, TNXB, HMCN2, MEP1A, SLC37A2, CHGB 
- Representative plotting workflows, stability summaries, and representative-seed display workflows are included in the testing scripts used during method development.

Current scope and limitations
StaBiCut is currently designed for candidate biomarker prioritization, not for exhaustive genome-wide high-throughput production use without additional engineering.
Current limitations include:
- the implementation is optimized for a candidate-gene workflow rather than a packaged high-performance genome-wide pipeline; 
- some assumptions are currently tailored to TCGA-style sample labeling; 
- the current version uses fixed component weights rather than a learned weighting scheme; 
- missing scoring components can propagate to the composite score in the current runner; 
- bootstrap stability evaluates local reproducibility of the selected cutoff, not global optimality across all resampled datasets. 
Users should therefore interpret StaBiCut as a transparent, auditable prioritization framework under explicit decision rules.

Data availability
Because the full study datasets may be large and/or governed by external usage policies, this repository is intended to contain:
- the analysis code, 
- minimal example inputs, 
- workflow documentation, 
- the output dictionary, 
- and representative lightweight example outputs where appropriate. 
Large intermediate files, full bootstrap result archives, and manuscript-linked release snapshots should be deposited in Zenodo and linked from the GitHub release.

Code availability
All StaBiCut source code used for the present implementation is available in this repository. Release-tagged snapshots used for manuscript submission should also be archived in Zenodo to provide a citable DOI-linked record.
For manuscript submission, please cite both:
1.	the GitHub repository URL, and 
2.	the DOI-linked archived release. 

Suggested GitHub vs Zenodo split
GitHub
- source scripts 
- README 
- LICENSE 
- CITATION.cff 
- sessionInfo.txt 
- workflow PNG 
- minimal examples 
- output dictionary 
- light representative outputs if desired 
Zenodo
- GitHub release snapshot 
- full results_seed_*.rds files 
- stability summary Excel files 
- representative-seed Excel files 
- other large manuscript-linked output archives 

Citation
If you use StaBiCut in your work, please cite:
- the associated manuscript, when available, 
- the GitHub repository, 
- and the archived Zenodo release corresponding to the version used. 
Example:
StaBiCut v2.0.0. GitHub repository and Zenodo archive.

Versioning
This repository follows semantic-style release tracking for manuscript-linked snapshots.
Suggested release for the present implementation:
- v2.0.0 
This label corresponds to the current StaBiCut implementation used for the reported study analyses and repository release preparation.

License
An explicit open-source license should be included before public release.
The current recommendation for this repository is:
- MIT License 

Contact
For questions, bug reports, or feature requests, please open a GitHub issue in this repository.



# Curated paper evidence

Results selected for the PRS-MIX manuscript
([statsleelab/PRSMIX_manuscript](https://github.com/statsleelab/PRSMIX_manuscript)).
Everything here is frozen and committed. The manuscript imports from this
folder only, through its explicit import step. It never reads working outputs
(`inst/scripts/results/`), scratch directories or cluster outputs.

One folder per study:

```text
paper/<study>/
  inputs/          frozen copies of the metric tables the assets are computed from
  *.tex            generated LaTeX table bodies (no caption; the manuscript adds it)
  tables.html      browser preview of the same tables (open locally; not imported)
  tables.tsv       every number in the tables, unrounded, with n per cell
  manifest.json    program, command, commits, settings, input and output sha256
```

Rules:

- Generate assets with a committed program, from committed inputs. Commit the
  program first, run it, then commit the outputs, so `manifest.json` records
  the exact program commit (`asset_commit`) and the commit of the inputs
  (`experiment_commit`).
- Never edit generated files by hand. Change the program or inputs and
  regenerate.
- Keep panels, genotypes, posterior draws and logs out of this folder.

## Studies

| Study | Program | Contents |
| --- | --- | --- |
| `null_sim_eur_asn_250kb` | `inst/scripts/null_sim/12_paper_tables.R` | chr22, EUR+ASN 70:30 GWAS, 250 kb blocks, 100 null and 100 non-null replicates, phi auto: null and non-null tables over the full weight sweep, paired differences against the PRS-CS seed rerun, and settings (LaTeX and HTML) |
| `sim250_sweep_250kb` | `inst/scripts/null_sim/18_sweep_paper_tables.R`; `figures/` by `19_sweep_paper_figures.R` (`footnote=no`) | chr22, 250 kb blocks, 100 replicates per setting, phi auto: accuracy (r^2, MSE), top-decile overlap against the seed rerun, multi-cutoff overlap and settings tables over 4 ancestry pairs x 3 GWAS compositions (the tables still include ASN+AFR, which the paper dropped on 2026-11-12; regenerate for the three EUR pairs before the manuscript imports them); `figures/` holds the EUR-pair top-10% overlap figure with its own manifest. Replaces the 20-replicate provisional placeholder, removed 2026-11-13. |
| `null_sim_ld_truth_250kb` | `figures/` by `inst/scripts/null_sim/21_ld_truth_figures.R` (`footnote=no`) | chr22, 250 kb blocks, block 75 (34.55-34.80 Mb): the LD the null GWAS z-scores were simulated with (discovery LD blended at the GWAS composition) against PRS-MIX's blended reference-panel LD and the majority (EUR) panel's LD, with both differences, at 90:10, 70:30 and 50:50, for EUR+ASN, EUR+AFR and EUR+SAS. No PRS-CS output and no replicates. `block_summary.tsv` has the mean \|difference\| values and the block's rank among all chr22 blocks. Inputs are the RedHawk extracts from `20_ld_truth_extract.R`. |
| `layout_pilot` | `inst/scripts/null_sim/22_pilot_paper_tables.R` | Supplementary tables from the prespecified layout pilot (`inst/scripts/null_sim/LAYOUT_PILOT.md`; EUR+AFR and EUR+ASN 50:50, chr22, replicates 1-20, chain arms 1-5) and the 100-replicate EUR+ASN 70:30 chain check: sampler chain length (S6a, S6b), reference vs discovery-built panels (S7), accuracy by block layout (S8a) and sampler-free LD diagnostics (S8b). 95% paired t intervals; LaTeX shows 250 kb, `tables.tsv` every layout. |
| `panel_validation` | `inst/scripts/null_sim/23_panel_validation_tables.R` | Genome-wide reference panels (chr1-22, EUR, ASN, AFR, AMR, SAS; 250 kb and PRS-CS EUR cut-point layouts): structural and matrix checks and agreement with PRS-CS's published 1000 Genomes panels (S1a), reference sample sizes and LD/MAF features (S1b), and PRS-MIX's script with one panel against unmodified PRS-CS (S1c). |

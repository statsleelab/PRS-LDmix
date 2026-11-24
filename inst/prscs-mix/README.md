# PRS-CS reference-panel mixing scripts

`PRScs_mix.py` and `v2parse_genet_mix.py` extend [PRS-CS](https://github.com/getian107/PRScs)
(Ge et al. 2019) to compute posterior SNP effects from a weighted blend of multiple
reference LD panels (e.g. EUR + ASN), instead of a single population panel. They are
not part of upstream PRS-CS and were not previously under version control anywhere.

To run, place these alongside an unmodified PRS-CS checkout so `mcmc_gtb.py` and
`gigrnd.py` are importable (`PRScs_mix.py` imports `v2parse_genet_mix`, `mcmc_gtb`,
`gigrnd`).

**License.** The two scripts are modified from PRS-CS, which is released under
the MIT License (copyright (c) 2018 Tian Ge); that notice is kept in
[`LICENSE-PRScs`](LICENSE-PRScs). The rest of PRS-MIX is released under the GNU
GPL, version 2 or later (see the repository's `LICENSE.md`).

`PRScs_mix.py`'s `main()` intersects every `--ref_dir` panel's SNP universe up front
before building the SNP list. Earlier versions built the SNP list from only the
first-listed panel, so a SNP missing from a smaller panel (e.g. ASN, which covers
fewer chr22 SNPs than EUR due to population-specific MAF filtering) was silently
dropped inside `parse_ldblk_mix`'s per-block LD computation instead of being excluded
from the SNP list -- leaving it with a posterior beta of exactly 0.0. Because the
sampler's shared noise-scale parameter in `mcmc_gtb.py` is drawn using the total SNP
count, those phantom zeroed SNPs also biased every other SNP's estimate.

`PRScs_mix.py` also checks `--mix_props` before running: one proportion per `--ref_dir`
panel, each non-negative, summing to 1 (within 1e-6). Earlier versions used the values as
given. `parse_ldblk_mix` pairs proportions with panels through `zip`, so a list shorter
than `--ref_dir` silently dropped the remaining panels, and proportions not summing to 1
gave a mixed LD matrix without a unit diagonal; neither raised an error.

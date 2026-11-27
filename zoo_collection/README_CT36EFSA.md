# CT36 EFSA: DEB models for food reduction on focal species

This folder holds the code and results behind the EFSA report
*EOI_EFSA_2022_01_CT 36 PLANTS, Revision 07_2026, TASK 1-4: DEB models for food
reduction on focal species* (version of 2026-11-23). The rest of the repository is a fork of
[add-my-pet/DEBtool_M](https://github.com/add-my-pet/DEBtool_M). Outside `zoo_collection/`,
the only change is a bug fix in `animal/get_tx.m` (`direction` of the `event_xp` events had
one element instead of two; the two shrews need it).

The report PDF is not part of the repository.

## Layout

```
zoo_collection/
  <Species>/                       9 species
    mydata_<sp>.m, pars_init_<sp>.m, predict_<sp>.m, run_<sp>.m   Add-my-Pet entry (unchanged)
    pars_init_<sp>_CT36EFSA.m      parameters + baseline f of each stage (f_base_*)
    predict_<sp>_CT36EFSA.m        predictions under a reduced f for one stage
  CT36EFSA_pipeline/               the only way to reproduce the report figures
    run_all_CT36EFSA.sh            everything: 9 species x 2 scenarios, then the figures
    run_octave_CT36EFSA.m          simulations of one species (GNU Octave)
    plot_CT36EFSA.R                figures of one species (R)
    event_CT36EFSA.m               age at a maturity event, NaN if not reached
    el_path_CT36EFSA.m             reserve and length along a piecewise-constant food path
    octave_shims/                  MATLAB functions missing or different in Octave
  CT36EFSA_results/                outputs of run_all_CT36EFSA.sh (figures + CSV)
  CT36EFSA_analysis/               Chapter 2 sensitivity analysis (R Markdown + CSV inputs)
```

## Scenarios

For each species and each life stage (nest or lactation or embryo, juvenile, adult), the
f of that stage is multiplied by a factor from 0.5 to 1.0 (step 0.01). The other stages
keep their baseline f.

- **acute** (report section 3.3): the reduced f lasts an n-day window that starts at the
  onset of the stage. The onsets are the baseline events (unrestricted f). The window is
  7 days for the birds, the shrews and *Podarcis*, 126 days for *Rana* and 336 days for
  *Ichthyosaura*.
- **full** (report Appendix C): the reduced f lasts the whole stage, which ends at the
  event reached at the reduced f. If that event is never reached, the restriction never ends.
  This scenario is selected with `par.ct36_full = 1` in the `predict_*_CT36EFSA.m` files.

The chronic outputs (`outputs_*.csv`: tp, Wwi, Ri, ...) hold f constant over the stage. They
are the same for both scenarios and are written once.

## Which file is which figure

| Report | File in `CT36EFSA_results/<species>/` |
|---|---|
| Fig 3.2 | `Alauda_arvensis_f_juvenile.png` |
| Figs 3.3-3.8, top | `<sp>_growth_combined.png` (Perdix, Motacilla, Alauda, Turdus, Sorex, Crocidura) |
| Fig 3.9, top | `Ichthyosaura_alpestris_growth_combined.png` (female row, male row) |
| Fig 3.10, top | `Rana_temporaria_growth_combined.png` + `Rana_temporaria_LdL_combined.png` |
| Fig 3.11, top | `Podarcis_muralis_combined_dynamics.png` (+ `_combined_fecundity.png`) |
| Figs 3.3-3.11, middle | `<sp>_growth_combined_relative.png` (Podarcis: `_combined_dynamics_relative.png`) |
| Figs 3.3-3.11, bottom | `<sp>_comparison_variables.png` |
| Figs C.1-C.6 | `<sp>_full_growth_combined.png` (Perdix, Motacilla, Alauda, Turdus, Sorex, Crocidura) |
| Fig C.7 | `Ichthyosaura_alpestris_full_growth_combined.png` |
| Fig C.8 | `Rana_temporaria_full_LdL_combined.png` |
| Fig C.9 | `Rana_temporaria_full_growth_combined.png` |
| Figs C.10, C.11 | `Podarcis_muralis_full_combined_dynamics.png` (SVL panels for C.10, weight-SVL panels for C.11) |
| Fig C.12 | `Podarcis_muralis_full_combined_fecundity.png` |
| Figs 2.2-2.17 | rendered by `CT36EFSA_analysis/dataAnalysis_CT36EFSA.Rmd` |

The CSV files behind the figures are in the same folders:
- `traj_<sp>_<stage>.csv` and `traj_<sp>_<stage>_full.csv`: growth curves (f_value, panel, x, y)
- `outputs_<sp>_<stage>.csv`: chronic outputs, one row per factor. An empty cell means that the
  event is not reached.
- `data_<sp>.csv`: Add-my-Pet data points drawn on the growth panels

`MANIFEST_CT36EFSA.sha256` holds the checksums of all result files.

## Reproduce

Requirements:
- GNU Octave 8.4. It was installed with `micromamba create -n octave -c conda-forge octave`.
  MATLAB is not used. The shims in `CT36EFSA_pipeline/octave_shims` are added automatically.
- R with `ggplot2`, `patchwork` and `viridis` for the figures. The Chapter 2 Rmd also needs
  `readr`, `dplyr`, `tidyr`, `gridExtra` and `broom`.

```bash
# from the repository root; about 30 min per scenario with 3 Octave processes
bash zoo_collection/CT36EFSA_pipeline/run_all_CT36EFSA.sh            # -> zoo_collection/CT36EFSA_results
OCTAVE=/path/to/octave JOBS=3 bash zoo_collection/CT36EFSA_pipeline/run_all_CT36EFSA.sh /tmp/check "acute"

# check a rerun against the committed results (CSV only: PNG bytes depend on the R version)
cd zoo_collection/CT36EFSA_results && sha256sum -c MANIFEST_CT36EFSA.sha256 2>/dev/null | grep -v ': OK' | grep csv

# Chapter 2 figures
cd zoo_collection/CT36EFSA_analysis && Rscript -e 'rmarkdown::render("dataAnalysis_CT36EFSA.Rmd")'
```

The Chapter 2 inputs (`ref_values - <sp>.csv`, `variation_f_food - <sp>.csv`) are the
0-variate predictions of the Add-my-Pet models for f between 1 and the lowest f the
model accepts (0.65 to 0.93, depending on the species). They were exported
by hand during Task 2; no script in this repository writes them.

## Version used in the report

The report of 2026-11-23 corresponds to the git tag `CT36EFSA-v11-20261123`.

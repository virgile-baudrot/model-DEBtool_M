DEBtool_M is Matlab code that can be used for application of Dynamic Energy Budget (DEB) theory
https://bio.vu.nl/thb/deb/deblab/

add-my-pet.github.io/debtool_m/docs/

## This fork: `zoo_collection`

Everything in this repository is a fork of
[add-my-pet/DEBtool_M](https://github.com/add-my-pet/DEBtool_M), except the
`zoo_collection/` folder and a bug fix in `animal/get_tx.m`. For DEBtool itself, see the
original repository and its documentation.

`zoo_collection/` holds the code, results and figures of the EFSA CT36 report on food
reduction in 9 focal species (4 birds, 2 shrews, 2 amphibians, 1 lizard):

- `<Species>/`: the Add-my-Pet entry of each species and its `_CT36EFSA` parameter and
  prediction files
- `CT36EFSA_pipeline/`: GNU Octave + R scripts that rerun all simulations and figures
- `CT36EFSA_results/`: figures and CSV outputs of the report (section 3.3 and Appendix C)
- `CT36EFSA_analysis/`: sensitivity analysis of the report's Chapter 2

The details are in [zoo_collection/README_CT36EFSA.md](zoo_collection/README_CT36EFSA.md):
which file is which report figure, and how to reproduce everything.

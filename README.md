[![Paper](https://img.shields.io/badge/paper-arXiv%3AXXXX.YYYYY-B31B1B.svg)](https://arxiv.org/abs/XXXX.YYYYY)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7271852.svg)](https://doi.org/10.5281/zenodo.7271852)

# Strain-induced superfluid transition for atoms on graphene

Sang Wook Kim, Mohamed Elsayed, Taras Lakoba, Juan Vanegas, Valeri Kotov, Adrian Del Maestro

[arXiv:XXXX.YYYYY](https://arxiv.org/abs/XXXX.YYYYY)

### Abstract
Bosonic atoms deposited on atomically thin substrates represent a playground for exotic quantum many-body physics due to the  highly-tunable, atomic-scale nature of the interaction potentials.  The ability to engineer strong interparticle interactions can lead to the emergence of complex collective atomic states of matter, not possible in the context of dilute Bose gases confined by optical lattices. While it is known that the first layer of adsorbed helium on graphene is permanently locked into a solid phase, we show by a combination of quantum Monte Carlo and mean-field techniques, that  simple isotropic (graphene) lattice expansion effectively unlocks a large variety of two-dimensional ordered commensurate, incommensurate, cluster atomic solid, and superfluid states for adsorbed atoms.  It is especially significant that an atomically thin superfluid phase of matter emerges under experimentally feasible strain values, with potentially supersolid phases in close proximity on the phase diagram.

### Description
This repository includes links, code, scripts, and data to generate the figures in a paper.

### Requirements
The data in this project was generated via path integral Monte Carlo (PIMC) simulation. You can find the source code for PIMC on https://github.com/DelMaestroGroup/pimc and raw data on Zenodo https://doi.org/10.5281/zenodo.7271852. Python notebooks and modules for analysis are in the [src](https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/tree/main/src) (See README.md in the directory), and data you need is in the [data](https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/tree/main/data) directory

* Dependency: See [./src/README.md](https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/tree/main/src).

You can also install a minimal environment via: `pip install -r requirements.txt` 

### Support

This work was supported by NASA grant number [80NSSC19M0143](https://www.usaspending.gov/award/ASST_NON_80NSSC19M0143_8000).

### Figures

#### Figure 01: Strain-tuning the mean field phase diagram
<img src="https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/blob/main/figures/mf_phase_diagrams.svg" width="400px">

#### Figure 02: Strain-dependent adsorption and model parameters
<img src="https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/blob/main/figures/VG_vs_z_BH.svg" width="400px">

#### Figure 03: Superfluid phase diagram for helium adsorbed on strained graphene
<img src="https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/blob/main/figures/fig3.svg" width="400px">

#### Figure 04: Details of the superfluid phase
<img src="https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/blob/main/figures/fig4.svg" width="400px">

#### Figure 05: Mean-field phase diagram in physical units driven from Hartree-Fock based model parameters
<img src="https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/blob/main/figures/mu_delta_phase_diagram_withS_high.svg" width="400px">

#### Figure 06: Example of finite size effect analysis
<img src="https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/blob/main/figures/figS3.svg" width="400px">

#### Figure 07: Averaged linear density along z-direction
<img src="https://github.com/DelMaestroGroup/papers-code-Superfluid4HeStrainGraphene/blob/main/figures/linden15-93.svg" width="400px">

This figures are relesed under [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/) and can be freely copied, redistributed and remixed.

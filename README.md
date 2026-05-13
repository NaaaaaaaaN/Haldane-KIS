# Haldane-KIS

This repository contains the data and code for the paper *The unseen complexity of inhibition kinetics: the Haldane KIS*, available in MATLAB.

## Data

The data is extracted from [1] and [2] for eight microalgae species. Each dataset contains two columns: one for light irradiance and one for growth rate. The data files are named as follow:
* Anning_PcHL = *Skeletonema costatum I<sub>high</sub>*,
* Anning_PcLL = *Skeletonema costatum I<sub>low</sub>*,
* Yang1 = *Isochrysis galbana*,
* Yang2 = *Dunaliella salina*,
* Yang3 = *Platymonus subcordiformis*,
* Yang4 = *Chlorococcum sp. FACHB-1556*,
* Yang5 = *Microcystis aeruginosa FACHB-905*,
* Yang6 = *Microcystis wesenbergii FACHB-1112*,
* Yang7 = *Scenedesmus obliquus FACHB-116*.

## Descriptions

Five descriptions of the Haldane model are considered: four existing in the literature and one newly proposed in the paper.

## MATLAB version

MATLAB version uses the function `fminsearch` in `Optimization Toolbox` to find the optimal parameters.

A single `main.m` script contains all necessary sub-functions to compute results and generate figures. To run the script `main.m`, the `data` folder also needs to be downloaded, and the data name needs to be specified on the second line of the script, e.g., `dataName='Anning_PcHL'`. The script provides all necessary results presented in the paper for the five descriptions. It also generates and saves in `fig` folder all figures illustrated in the Appendix.

## References

[1] T. Anning, H.L. MacIntyre, S.M. Pratt, P.J. Sammes, S. Gibb, and R.J. Geider. Photoacclimation in the marine diatom skeletonema costatum. Limnology and Oceanography, 45(8):1807–1817, 2000.

[2] X.L. Yang, L.H. Liu, Z.K. Yin, X.Y Wang, S.B Wang, and Z.P. Ye. Quantifying photosynthetic performance of phytoplankton based on photosynthesis–irradiance response models. Environmental sciences Europe, 32(1), 2020.

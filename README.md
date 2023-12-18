# Plasmid-mediated phenotypic noise leads to transient antibiotic resistance in bacterias

 Scripts and data necessary to produce figures and movies from

**Plasmid-mediated phenotypic noise leads to transient antibiotic resistance in bacterias**\
JCR Hernandez-Beltran, J Rodriguez-Beltran, B Aguilar Luviano, J Velez-Santiago, O Mondragon-Palomino, RC MacLean, A Fuentes-Hernandez, A San Millan, R Peña-Miller.\


## Overview

This supplementary material accompanies our study on the role of multicopy plasmids in transient antibiotic resistance in bacteria. The focus is on a computational model that integrates both single-cell and population-level dynamics to explore how variability in plasmid copy number (PCN) impacts bacterial survival under antibiotic pressure. Our model highlights the phenomenon of phenotypic noise resulting from PCN variability, offering insights into the survival strategies of bacterial populations. Through simulations, we demonstrate the impact of PCN diversity on bacterial adaptability and resistance, underlining the role of multicopy plasmids in these processes.

## Figures

**Figure 1:** Stochastic plasmid dynamics yield heterogeneous populations.\
Jupyter Notebook: [code/Fig_1-A_jl.ipynb](Fig_1-A_jl.ipynb)
Jupyter Notebook: [code/Fig_1-BCDEF_jl.ipynb](Fig_1-BCDEF_jl.ipynb)

**Figure 2:** Experimental model system.\
Jupyter Notebook: (...)

**Figure 3:** Effect of antibiotics on GFP and survival under fluctuating selection.\
Jupyter Notebook: (...)

**Figure 4:** Single-cell analysis of a semi-lethal pulse.\
Jupyter Notebook: (...)

**Figure 5:** Microscopy montage of a microfluidics semi-lethal pulse.\
Jupyter Notebook: (...)

## Supplementary Information: Agent-based model

The agent-based model simulates the dynamics of bacterial populations with multicopy plasmids under various antibiotic pressures. It incorporates key biological processes like plasmid replication, bacterial growth, and antibiotic-induced mortality. Through computer simulations, we explore how plasmid copy number variability influences the survival and adaptability of bacterial populations, providing a computational framework to better understand the role of plasmids in antibiotic resistance.

PDF: [Agent-Based Model](SI_pBGT_ABM.pdf)
Jupyter Notebook: (...)

## Supplementary Figures

PDF: [Supplementary Figures](SI_pBGT.pdf)\
Jupyter Notebook: (...)

## Supplementary Movies

Jupyter Notebooks used to generate movies are coded in
[Napari](https://napari.org/) and can be found in
[pBGT/code/Supplementary/](/code/Supplementary)

**Movie S1:** MG/pBGT exposed to an antibiotic ramp.\
Cell viewer:
[SuppMovie1.ipynb](/code/Supplementary/SuppMovie1.ipynb)\
AVI:
[SuppMovie1.avi](/movies/SuppMovie1.avi)

**Movie S2:** MG/pBGT exposed to a semi-lethal pulse of AMP.\
Cell viewer:
[SuppMovie2.ipynb](/code/Supplementary/SuppMovie2.ipynb)\
AVI:
[SuppMovie2.avi](/movies/SuppMovie2.avi)


## Authors
[\@Systems Biology Lab, CCG-UNAM](https://github.com/ccg-esb-lab)



## License

[MIT](https://choosealicense.com/licenses/mit/)

This project is licensed under the MIT License - see the
[license.txt](/ccg-esb-lab/pBGT/blob/main/license.txt) file for details.


# Plasmid-mediated phenotypic noise leads to transient antibiotic resistance in bacterias

 Scripts and data used to produce figures and movies from

**Plasmid-mediated phenotypic noise leads to transient antibiotic resistance in bacterias**\
JCR Hernandez-Beltran, J Rodriguez-Beltran, B Aguilar Luviano, J Velez-Santiago, O Mondragon-Palomino, RC MacLean, A Fuentes-Hernandez, A San Millan, R Peña-Miller.


## Overview

This is the repositofy of our study on multicopy plasmid heterogeneity and its role in transient antibiotic resistance in bacteria. It containes a series of scripts for data analysis at both single-cell and population levels. Central to our study is a computational model that evaluates into the interplay between individual cell dynamics and collective bacterial behavior. This model underscores the significance of plasmid copy number variability and its resultant phenotypic noise, offering deep insights into bacterial survival mechanisms under antibiotic stress and the adaptive capabilities conferred by multicopy plasmids.

## Figures

**Figure 1:** Stochastic plasmid dynamics yield heterogeneous populations.\
Jupyter Notebook: [Fig_1-A_jl.ipynb](code/Fig_1-A_jl.ipynb)\
Jupyter Notebook: [Fig_1-BCDEF_jl.ipynb](code/Fig_1-BCDEF_jl.ipynb)

**Figure 2:** Experimental model system.\
Jupyter Notebook: [py-Fig_2CDE.ipynb](code/py-Fig_2CDE.ipynb)

**Figure 3:** Effect of antibiotics on GFP and survival under fluctuating selection.\
Jupyter Notebook: [py-Fig_3.ipynb](code/py-Fig_3.ipynb)

**Figure 4:** Single-cell analysis of a semi-lethal pulse.\
Jupyter Notebook: [py-Fig_4BCDEF.ipynb](code/py-Fig_4BCDEF.ipynb)

**Figure 5:** Microscopy montage of a microfluidics semi-lethal pulse.\
Jupyter Notebook: [py-Fig_5CD.ipynb](code/py-Fig_5CD.ipynb)

## Supplementary Figures

PDF: [Supplementary Figures](SI_pBGT.pdf)\
Jupyter Notebook: [py-Fig_S1.ipynb](code/Supplementery/SupplementaryFigures/py-Fig_S1.ipynb)\
Jupyter Notebook: [py-Fig_S2.ipynb](code/Supplementery/SupplementaryFigures/py-Fig_S2.ipynb)\
Jupyter Notebook: [py-Fig_S3.ipynb](code/Supplementery/SupplementaryFigures/py-Fig_S3.ipynb)\
Jupyter Notebook: [py-Fig_S4.ipynb](code/Supplementery/SupplementaryFigures/py-Fig_S4.ipynb)\
Jupyter Notebook: [py-Fig_S5.ipynb](code/Supplementery/SupplementaryFigures/py-Fig_S5.ipynb)\
Jupyter Notebook: [py-Fig_S6-S7-S8-S9-S10-S11-S12.ipynb](code/Supplementery/SupplementaryFigures/py-Fig_S6-S7-S8-S9-S10-S11-S12.ipynb)\
Jupyter Notebook: [py-Fig_S13.ipynb](code/Supplementery/SupplementaryFigures/py-Fig_S13.ipynb)

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


## Supplementary Information: Agent-based model

The agent-based model simulates the dynamics of bacterial populations with multicopy plasmids under various antibiotic pressures. It incorporates key biological processes like plasmid replication, bacterial growth, and antibiotic-induced mortality. We use computer simulations implemented in Julia to explore how plasmid copy number heterogeneity influences the survival of bacterial populations, providing a computational framework to better understand the role of plasmids in the evolution of antibiotic resistance.

PDF: [Agent-Based Model](SI_pBGT_ABM.pdf)
Jupyter Notebook: [SupplementaryInformation.ipynb](code/Supplementary/SupplementaryInformation/SupplementaryInformation.ipynb)

## Authors
[\@Systems Biology Lab, CCG-UNAM](https://github.com/ccg-esb-lab)



## License

[MIT](https://choosealicense.com/licenses/mit/)

This project is licensed under the MIT License - see the
[license.txt](/ccg-esb-lab/pBGT/blob/main/license.txt) file for details.


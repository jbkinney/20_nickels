## Computational supplement for Skalenko et al. (2021)

This repository contains analysis scripts and processed data for, 

**Skalenko et al. (2021) "Promoter sequence determinants and structural basis of primer dependent transcription initiation in *Escherichia coli*".  [bioRxiv doi:10.1101/2021.04.06.438613](https://www.biorxiv.org/content/10.1101/2021.04.06.438613), Proc. Natl. Acad. Sci. USA. *In press*.**

To reproduce the computational analyses in this manuscript, execute the following iPython notebooks in this order.

- 1_efficiency.ipynb
- 2_mixture_modeling.ipynb
- 3_compute_invivo_logos.ipynb
- 4_compute_invitro_logos.ipynb
- 5_compute_chromosomal_logo.ipynb

Using pre-processed data in `data/`, these scripts will perform analyses and write the results to CSV files in `csv_results/` and `csv_logos/`. To then create select figure panels, execute the following iPython notebooks.

- fig_3A_oh_pct_table.ipynb
- fig_4_invivo_invitro_logos.ipynb
- fig_6A_chromosomal_logo.ipynb
- fig_S1A_mixture.ipynb
- fig_S5_invivo_logos_by_primer.ipynb
- fig_S6_low_hi_invitro_logos.ipynb
- fig_S8_profiles.ipynb

The resulting graphics will be written to `figures/`. Note that all output directories are already populated with the output of these analyses. 

Please address technical questions about this repository and its contents to [Justin B. Kinney](mailto:jkinney@cshl.edu). More general scientific correspondence about this work should be sent to [Bryce Nickels](mailto:bnickels@waksman.rutgers.edu).

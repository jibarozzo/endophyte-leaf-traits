# Evaluating Endophyte-Rich Leaves and Leaf Functional Traits’ Role in Protecting Tropical Trees Against a Generalist Herbivore and a Pathogen

In this study, we investigate the role of FEF abundance, diversity and community composition in the modulation of leaf functional traits and plant's response to herbivory and pathogen damage. The dataset includes information about anthocyanin content, leaf thickness, toughness, and leaf mass per area for different tropical woody species and the foliar endophytic fungi found inside their leaf tissue. 

This repository contains the R code and notebooks used for bioinformatic pipeline and statistical analyses used in this project. We include a first manuscript (*Aponte_Bolivar_2023_Aim3_manuscript.docx*) and the main notebook used for analyses (*Aim3_leaf_traits_statistical_analyses_v3.qmd file*). Supplementary files are also included in the repository.

>[!IMPORTANT]
>As of 2024-09-09 changes to the manuscript content will be completed in the .docx file. The `Aponte_Bolivar_2023_Aim3_manuscript.qmd` and `Aim3_leaf_traits_statistical_analyses_v3.qmd` file are almost 100% updated in par with .docx but some differences may persist.

## Table of Contents

1. [Project Structure](#project-structure)
2. [Contributing](#contributing)

# Project Structure
`manuscript`: Reproducible manuscript output.

`field_data/`: Contains the raw data collect in field experiments.

`clean_data/`: Contains the cleaned microbial data frames used for taxonomic, phylogenetic, and statistical analyses. Contains the R objects generated during the analysis. Used for accelerating the project workflow after completion of analyses.

`post_stat_analysis_data/`: Contains the OTU data obtained after `multipatt` analysis and selection cut-offs.

`filters/`: Contains the filter joins for sub-setting data.

`tables/`: Output directory for generated tables.

`figures/`: Output directory for generated figures.

# Contributing
Feel free to contribute by opening issues, providing feedback, or submitting pull requests.


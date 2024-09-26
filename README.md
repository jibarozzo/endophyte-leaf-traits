# Evaluating Endophyte-Rich Leaves and Leaf Functional Traits’ Role in Protecting Tropical Trees Against a Generalist Herbivore and a Pathogen

In this study, we investigate the role of FEF abundance, diversity and community composition in the modulation of leaf functional traits and plant's response to herbivory and pathogen damage. The dataset includes information about anthocyanin content, leaf thickness, toughness, and leaf mass per area for different tropical woody species and the foliar endophytic fungi found inside their leaf tissue.

This repository contains the R scripts, notebooks and sequences used for bioinformatic pipeline and statistical analyses used in this project. The DNA sequences for all the samples can be obtained in the NCBI BioProject PRJNA1162076. We include a first manuscript (*Aponte_Bolivar_2023_Aim3_manuscript.docx*) and the main notebook used for analyses (*Aim3_leaf_traits_statistical_analyses_v3.qmd*). Supplementary files are also included in the repository.

*\[!IMPORTANT\] As of 2024-09-09 changes to the manuscript content will be completed in the .docx file. The `Aponte_Bolivar_2023_Aim3_manuscript.qmd` and `Aim3_leaf_traits_statistical_analyses_v3.qmd` file are almost 100% updated in par with `.docx` but some differences may persist.*

## Table of Contents

1.  [Project Structure](#project-structure)
2.  [Contributing](#contributing)
3.  [Reproducing analyses](#reproducing-analyses)

# Project Structure {#project-structure}

Main directories and files in the project:

```         
├── Aim3_endophytes_leaf_traits.Rproj
├── Aim3_leaf_traits_statistical_analyses.qmd
├── Aponte_Bolivar_aim3_manuscript.qmd
├── clean_data
├── field_data
├── figures
├── filters
├── functions
├── manuscript
├── mock_ouputs
├── ncbi
├── post_stat_analyses_data
├── references.yaml
├── renv
├── renv.lock
├── scripts
├── tables
└── unite_report2BSDHQZP.csv
```

`manuscript`: Reproducible manuscript output. `field_data/`: Contains the raw data collect in field experiments. `clean_data/`: Contains the cleaned microbial data frames used for taxonomic, phylogenetic, and statistical analyses. Contains the R objects generated during the analysis. Used for accelerating the project workflow after completion of analyses. `post_stat_analysis_data/`: Contains the OTU data obtained after `multipatt` analysis and selection cut-offs. `filters/`: Contains the filter joins for sub-setting data. `tables/`: Output directory for generated tables. `figures/`: Output directory for generated figures. `scripts/`: Contains the R scripts used for the bioinformatic pipeline and statistical analyses. `functions/`: Contains the R functions used in the scripts. `ncbi/`: Contains the sequences used for the bioinformatic pipeline. *See NCBI BioProject PRJNA1162076 for more information.* `mock_outputs/`: Contains the mock outputs for the pipeline. `unite_report2BSDHQZP.csv`: The UNITE database report used for the taxonomic assignment.

*File paths and repo structure might need to be updated to you local machine*

## Reproducing analyses {#reproducing-analyses}

1.  Copy repository to your local machine
2.  Download the sequences from the NCBI BioProject PRJNA1162076
3.  Put the sequences in the `ncbi/` directory
4.  Run the `V-Search_Post_Sequencing.md` pipeline in the `scripts/` directory
5.  Run the `Aim3_leaf_traits_statistical_analyses.qmd` notebook.

-   This notebook contains more bioinformatic cleaning steps and the statistical analyses. It is ordered by analyses types.

6.  Run the `Aponte_Bolivar_aim3_manuscript.qmd` notebook.

-   This notebook contains the manuscript content and the figures and tables generation. It contains refined analyses and results.

# Contributing {#contributing .contributing}

Feel free to contribute by opening issues, providing feedback, or submitting pull requests.

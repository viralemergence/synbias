# Coinfection Dynamics and Community Surveillance

This repository contains the data and code for the manuscript by Pilowsky *et al.* on how coinfection dynamics influence observations of parasite and pathogen communities. The study employs both theoretical simulations and empirical data analysis to explore the interactions between multiple pathogens within host communities.

## File Structure

```
.
├── _quarto.yml
├── Data
│   ├── Empirical
│   │   ├── Parasite_data.xlsx
│   │   ├── Trematode Infections.csv
│   │   └── Trematode surveys at Carpinteria Salt Marsh, California USA, Feb 2012 to Jan 2014.xml
│   ├── epidemic_duration_seir.csv
│   ├── epidemic_duration_si.csv
│   ├── interaction_scores.csv
│   ├── round7_seir_output.csv
│   ├── round7_si_output.csv
│   ├── seir_cumulative_detection.csv
│   ├── simulation_round1.csv
│   ├── simulation_round2.csv
│   ├── simulation_round3.csv
│   ├── simulation_round4.csv
│   ├── simulation_round5
│   ├── simulation_round6
│   └── simulation_round7_input.csv
├── Figures
│   ├── conceptual_figure.pdf
│   ├── conceptual_figure.png
│   ├── deprecated
│   ├── Fig2a_seir_species.png
│   ├── Fig2a_si_species.png
│   ├── Fig2b_si_species.png
│   ├── Fig3_shapley.png
│   ├── Fig4_empirical_summary.png
│   ├── Fig4a_snail.png
│   ├── Fig4b_daphnia.png
│   ├── Fig5_species_interactions_curves.png
│   ├── interaction_accumulation_curve.png
│   ├── shapley_plots_seir.png
│   ├── shapley_plots_si.png
│   ├── supp_fig1a_errors.png
│   └── supp_fig1b_errors.png
├── README.md
├── Scripts
│   ├── _publish.yml
│   ├── deprecated
│   ├── empirical_analysis.qmd
│   ├── Figure scripts
│   │   ├── conceptual_figure.jl
│   │   ├── figure4.R
│   │   └── shapley_plots.R
│   ├── references.bib
│   └── theoretical_simulations.qmd
└── synbias.Rproj
```

## Description of Contents

- **Data/**: Contains all datasets used in the study, including empirical data from parasite surveys and simulation outputs from theoretical models. Round 7 is the most recent round of simulations; earlier rounds are included for transparency about our scientific process.
    - **Empirical/**: Subdirectory containing empirical datasets used for analysis in the manuscript.
- **Figures/**: Contains all figures generated for the manuscript, including main figures and supplementary materials. The "deprecated" folder includes older versions of figures from previous rounds of simulations.
- **Scripts/**: Contains all code used for data analysis and figure generation. The "Figure scripts" subdirectory includes scripts specifically for generating figures. The "deprecated" folder contains older versions of scripts from previous rounds of simulations.
  - **empirical_analysis.qmd**: Quarto document for analyzing empirical data. Written entirely in R; can be rendered to a HTML or PDF report.
  - **theoretical_simulations.qmd**: Quarto document for running and analyzing theoretical simulations. Written in a mix of Julia and R, and therefore cannot be rendered to a report, but can be run chunk by chunk on a machine with both Julia and R installed.

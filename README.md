# Introduction

This folder contains the code used in the Workplace Norms and Electric Vehicle Charging paper as well as the anonymized version of the data upon which the paper is built. It is designed for easy replication of all figures and tables.

# Code Instructions and Information

All figures and tables are generated from the single R file saved in this folder, "package.R". Running the full file sequentially will automatically generate and save all figures and tables that are derived from data (including all but Figure 1, which is text-based).

## Dependencies

The following R packages are used in this analysis:
1. "rddensity"
2. "ggplot2"
3. "rdd"
4. "car"
5. "dplyr"
6. "lubridate"
7. "dummies"
8. "lpdensity"
9. "extrafont" (optional for replication)

## Navigation

Once the code has executed, each table and figure from the paper can be found it its own variable.

The naming conventions for these variables are designed to be intuitive and searchable, with the following naming conventions:

For figures -- "figure_#"

For tables -- "table_#"

In cases where a figure has multiple panes, the format is as follow: "figure_{figure #}{pane letter}". For example, Figure 2A appears as "figure_2a".

Figures may be generated non-sequentially in some cases and are organized by the raw data used to generate each rather than by the order in which they appear in the paper. 





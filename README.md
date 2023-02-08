This repository contains public code for the manuscript:

Zou, Heng-Xing, and Volker HW Rudolf. "Priority effects determine how dispersal affects biodiversity in seasonal metacommunities." bioRXiv (2022). [https://doi.org/10.1101/2022.02.03.479022](https://doi.org/10.1101/2022.02.03.479022)

This manuscript has been accepted by *The American Naturalist*.

For any questions, please [contact Heng-Xing Zou](hengxingzou@rice.edu)

# Abstract

The arrival order of species frequently determines the outcome of their interactions. This phenomenon, called the priority effect, is ubiquitous in nature and determines local community structure, but we know surprisingly little about how it influences biodiversity across different spatial scales. Here, we use a seasonal metacommunity model to show that biodiversity patterns and the homogenizing effect of high dispersal depend on the specific mechanisms underlying priority effects. When priority effects are only driven by positive frequency dependence, dispersal-diversity relationships are sensitive to initial conditions but generally show a hump-shaped relationship: biodiversity declines when dispersal rates become high and allow the dominant competitor to exclude other species across patches. When spatiotemporal variation in phenological differences alters species’ interaction strengths (trait-dependent priority effects), local, regional, and temporal diversity are surprisingly insensitive to variation in dispersal, regardless of the initial numeric advantage. Thus, trait-dependent priority effects can strongly reduce the effect of dispersal on biodiversity, preventing the homogenization of metacommunities. Our results suggest an alternative mechanism that maintains local and regional diversity without environmental heterogeneity, highlighting that accounting for the mechanisms underlying priority effects is fundamental to understanding patterns of biodiversity.

# Code Files

All simulations were run in R version 4.2.1.

- `Functions_SpatialBH.R` contains all functions used for contructing the model and running simulations. 
- `FiguresReadMe.Rmd` contains all code for generating figures in the main text and the appendix. Running any portion of this file will automatically run `Functins_SpatialBH.R`. Note that running all simulations takes a long time; we therefore do not recommend knitting this file. Instead, we have provided a `.RData` file that contains all generated data.

# Data

- `SavedMetacoms` contains initial metacommunities in `.csv` files. You can also generate metacommunities in `FiguresReadMe.Rmd`, but resulting figures will likely be different because of the random generation process. 
- `Dec17_5spp.RData` contains generated data from all simulations in the main text. This is a pre-packaged R dataset that includes all variables in the environment. To load the dataset, simply double click while `Functions_SpatialBH.R` and `FiguresReadMe.Rmd` are opened in RStudio. Loading `Dec17_5spp.RData` will load the following key datasets:
  * `lowdisp_1`, `lowdisp_2`, `highdisp_1`, `highdisp_2`: Population dynamics of all patches and all repetitions. `low` or `high` denotes low or high dispersal rates. `_1` and `_2` denote scenarios (1. Equal Initials; 2. First Colonizer).
  * `num_Sce1`, `num_Sce2`, `tra_Sce1`, `tra_Sce2`: Alpha, beta and gamma diversity under a range of dispersal rates of all repetitions. `num` indicates the absence of trait-dependent priority effects, and `tra` indicates the presence. `Sce1` and `Sce2` refer to scenarios (1. Equal Initials; 2. First Colonizer).
  * `num_lowvar_1`, `tra_lowvar_1`, `num_lowvar_2`, `tra_lowvar_2`, and similarly named data frames containing `midvar` or `highvar` instead of `lowvar`: Alpha, beta and gamma diversity under a range of dispersal rates, three different variations of emergence times, and of all repetitions. `num` indicates the absence of trait-dependent priority effects, and `tra` indicates the presence. `lowvar`, `midvar`, and `highvar` denote the magnitude of variations. `_1` and `_2` denote scenarios (1. Equal Initials; 2. First Colonizer)
  * Loading `Dec17_5spp.RData` will also load many auxiliary data files, which are automatically generated when plotting and do not require re-running time-consuming simulations.

# Reproducing Figures

- Unarchive the provided `SavedMetacoms.zip`.
- Make sure all files, including the `SavedMetacoms` folder, are under the same directory.
- Open `FiguresReadMe.Rmd`, then load `Dec17_5spp.RData`.
- Follow specific instructions in `FiguresReadMe.Rmd` to generate figures in the main text, or change model parameters to produce figures in the appendix.

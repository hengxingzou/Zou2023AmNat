This repository contains public code for the manuscript:

Zou, Heng-Xing, and Volker HW Rudolf. "Priority effects determine how dispersal affects biodiversity in seasonal metacommunities." bioRXiv (2022). [https://doi.org/10.1101/2022.02.03.479022](https://doi.org/10.1101/2022.02.03.479022)

For any questions, please [contact Heng-Xing Zou](hengxingzou@rice.edu)

# Abstract

The arrival order of species frequently determines the outcome of their interactions. This phenomenon, called the priority effect, is ubiquitous in nature and determines local community structure, but we know surprisingly little about how it influences biodiversity across different spatial scales. Here, we use a seasonal metacommunity model to show that biodiversity patterns and the homogenizing effect of high dispersal depend on the specific mechanisms underlying priority effects. When priority effects are only driven by positive frequency dependence, dispersal-diversity relationships are sensitive to initial conditions but generally show a hump-shaped relationship: biodiversity declines when dispersal rates become high and allow the dominant competitor to exclude other species across patches. When spatiotemporal variation in phenological differences alters species’ interaction strengths (trait-mediated priority effects), local, regional, and temporal diversity are surprisingly insensitive to variation in dispersal, regardless of the initial numeric advantage. Thus, trait-mediated priority effects can strongly reduce the effect of dispersal on biodiversity, preventing the homogenization of metacommunities. Our results suggest an alternative mechanism that maintains local and regional diversity without environmental heterogeneity, highlighting that accounting for the mechanisms underlying priority effects is fundamental to understanding patterns of biodiversity.

# Code Files

- `Functions_SpatialBH.R` contains all functions used for contructing the model and running simulations.
- `FiguresReadMe.Rmd` contains all code for generating figures in the main text and the appendix. Running any portion of this file will automatically run `Functins_SpatialBH.R`. Note that running all simulations takes a long time; we therefore do not recommend knitting this file. Instead, we have provided a `.RData` file that contains all generated data.

# Data

- `SavedMetacoms` contains initial metacommunities in `.csv` files. You can also generate metacommunities in `FiguresReadMe.Rmd`, but resulting figures will likely be different because of the random generation process. 
- `Dec17_5spp.RData` contains generated data from all simulations in the main text. 

# Reproducing Figures

- Unarchive the provided `SavedMetacoms.zip`.
- Make sure all files, including the `SavedMetacoms` folder, are under the same directory.
- Open `FiguresReadMe.Rmd`, then load `Dec17_5spp.RData`.
- Follow specific instructions in `FiguresReadMe.Rmd` to generate figures in the main text, or change model parameters to produce figures in the appendix.

# Analysis of Hydrocarbonoclastic Microbial Communities Subjected to Different Oxygenation Conditions

## Abstract
This internship report presents the analysis of hydrocarbonoclastic microbial communities subjected to various 
oxygenation conditions, as part of research conducted at IPREM. The main objective is to develop and apply 
biostatistical and bioinformatics tools to analyze and characterize these microbial communities, 
focusing on sulfate-reducing bacteria (SRB) and denitrifying bacteria (DB). 
Three distinct microbial communities were incubated in bioreactors under conditions of permanent anoxic, 
permanent oxic, and anoxic/oxic oscillations. Sequencing data, processed with Qiime2 and statistically analyzed with 
Python, show that oxygenation conditions significantly influence community structure. DBs adapt better to oscillating 
conditions, while the abundance of SRBs decreases after oxygenation. The expression of dsrB genes and the relative 
abundance of ASVs indicate that certain bacteria, such as Sulfurimonas and Marinobacter, play crucial roles in these 
environments. These results provide important insights into the dynamics of microbial communities and their role in 
the bioremediation of contaminated marine sediments.

## Keywords
Microbial communities, biostatistics, bioinformatics, hydrocarbons

## Files
- `README.md`: this file
- `Rapport de stage.pdf`: the internship report (in french)
- `requirements.txt`: the list of required libraries
- `fig_2.a.py`: This script preprocesses data, calculates the Bray-Curtis distance between data points, performs Principal Coordinate Analysis (PCoA) on the initial conditions, and generates 3D plots and variance explained plots, saving them in various formats.
- `fig_2b.py`: This script preprocesses data, calculates the Bray-Curtis distance, performs PCoA, generates 3D plots and variance explained plots, saves them in various formats, and creates a rotating 3D plot gif.
- `fig_3.py`: This script contains functions for calculating the observed co-occurrence of species in a given presence-absence matrix. It also includes functions for loading data, calculating positions of species co-occurrence, and plotting the results.
- `fig_4.py`: This script contains functions for analyzing and visualizing the relative abundance of different bacterial species over time under different conditions. It also includes functions for performing statistical tests on the data.
- `data/explorer.py`: This script contains functions for loading and preprocessing the data. It combines several columns into a single 'Taxonomy' column and creates a dictionary of names for each unique taxonomy.

## Usage
To use these scripts, you need to have Python installed along with the necessary libraries,
described in the [requirements](requirements.txt).

First, load the data using the true_loader function from `data/explorer.py`.
This will return a pandas DataFrame with the preprocessed data and a dictionary of names for each unique taxonomy.
Then, you can use the functions in `fig_2_a.py`, `fig_2_b.py`, `fig_3.py` and `fig_4.py` to analyze and visualize 
the data.

I didn't included the data I used for this analysis for privacy reasons. But you can check the final results in the
internship report and the [final](final) folder.
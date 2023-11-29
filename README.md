# cadth2023-ROMA
CADTH Symposium Course: Using RWE for Research-Oriented Market Access for High-Cost Therapies

### Are you new to R?
R is one of the most popular programing languages in use, and has a massive community and open source custom packages. 
The official R software environment is a free open source environment with a rudamentary command line interface.
There is a rather ubiquitous third-party graphical user interface called 'R Studio' that significantly improves R's core functionality, a better command script console, quick-reference scripting, and (it must be said) aesthetics.
While you only require the R environment, we recommend you download R Studio as well.

Downloads:
R: https://cran.r-project.org/

R Studio: https://www.rstudio.com/products/rstudio/download/

Downloading help: https://rstudio-education.github.io/hopr/starting.html

While we will keep this directory available in the future, we recommend saving the files in this github directory in a folder that will also be your working directory. This is not required to run the script but will keep you organized with setting a working directory for your R workspace, which is necessary to run the script.
Working Directory = A file location you can set in R (command = setwd("File Path") where R will save your data (save()) and workspace (save.image())

Copy and save the full script of each file in this directory into its own R script.

Core files available:
### CADTH2023 ToyMarkov
The script for creating the base case Markov Model, with probabilistic sensitivity analysis and key results

In lines 8-15, we list a couple add-on packages we recommend you install in R Studio. To do so, all you need to do is run each of the following lines in a blank R script:

install.packages("dplyr")

install.packages("ggplot2")

'dplyr' is part of the 'tidyverse' package, and includes many custom functions that will make data manpulation easier

'ggplot2' is a visualization package that gives you a more straight forward set of graphing options to make pretty plots

Once installed, you only need to tell R to use access these add ons (library()) for you to have access to each package. The current script in this file assumes you have already installed each of these packages and only uses library().

### CADTH2023 ToyMarkov_Deterministic
A simplified version of the base case model that only generates deterministic values. Not worth much analytically, but potentially helpful if you're new to HTA in R

### CADTH2023 VOI
The script for performing Value of Information analysis on the results of the ToyMarkov model.
You must have the CADTH2023ToyMarkov R Workspace saved on the working directory to successfully run this script, or run the ToyMarkov.R file first, followed by this file (ignoring the load() command on line 6). 
We recommend you do not run the entire script at once since it will take quite some time to run (current clock = 1.4 hours).

[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/rAntonioh/AAAGs_2018/master)

# AGAR 2018 University of Buffalo
#### Author: R. Antonio Herrera
##### Contact: dr@ranton.io
##### Code created for this workshop is freely available to use/reuse
------------------
### Molecular genetic and genomic analysis of DNA mutations

##### Workshop:
- This is an interactive and guided workshop in three parts.
  1. Introduction

  2. Part 1. Source of variation and innovation.

  3. Part 2. Molecular genetic and genomic modeling of mutations.

  4. Part 3. Gene Ontology Enrichment Analysis

##### Learning Objectives:
- After attending this workshop, you will be able to:
  1. Explain how mutations are sources of variation in evolution & development.

  2. Create queries to NCBI for gene specific data.

  3. Simulate gene regulatory dynamics.

  4. Model molecular interactions.

  5. Perform a basic GO analysis.

##### Prerequisites:
1. Create a free account on http://www.pythonanywhere.com

2. Under Accounts, click the "teacher" tab and add "**rah**"" as your teacher.

3. ***NB!*** This code is written for Python3.6, compatibility has not been tested.

##### Background:
- Multicellular organisms are a form of complex life which develop from a single cell zygote to adult, under the control of the same DNA genome. Basic cellular processes include:
  1. Proliferation

  2. Differentiation into types

- The evolution of development studies how genomic mutations affect the changes in cellular complexity under stabilizing selection. For example, how do organisms evolve:
  1. Innovative and/or adaptive traits

  2. Increased mutational variation

  3. Gene network robustness against perturbation

##### Setup:
- To start to characterize how different mutations in genes affect the evolution of robustness in genomes of multicellular organisms we need to get things ready:
- Navigate to your Dashboard on http://www.pythonanywhere.com
  1. Open a bash console

  2. Make sure you are in your /home/username directory

  3. View file contents:

    <code>ls -lh</code>

  4. If you have both install.sh and a BioNetGen tarball run the command:

    <code>./install.sh</code>

  5. To test the installation go to your Files tab

  6. Open the epas1_biopython.py in the editor

  7. Read the code and run **once** in the console

### Part 1. Source of variation and innovation:
- Slides - How mutation leads to genetic diversity.
- Interactive gene regulatory network dynamics in python.

#### Instructions
- GRN modeling
  1. Navigate to https://mybinder.org/v2/gh/rAntonioh/AAAGs_2018/master

  2. Wait for the server to build.

  3. Enter Session6-GeneticAnalysis

  4. Start the GRN Model for AGAR-2018.ipynb Jupyter notebook by clicking it.

  5. After running the notebook, perform the Tasks.

### Part 2. Molecular genetic and genomic modeling of mutations:
- Slides - From structure to function and phenotype.
- Interactive rule-based modeling of molecular species and regulatory dynamics.

#### Instructions
- PythonAnywhere
  1. In chrome, open 2 tabs, one for a bash console and another to your files tab
  2. In the bash console, edit egfr_pysb.py using nano and enter your PythonAnywhere username
  3. In the bash console run the following:

    <code>python egfr_pysb.py</code>

  4. When complete, download to view egfr_observables.png from your Files tab
  5. Edit the matplotlib code to label the observables and give the figure a title.

### Part 3. Gene Ontology Enrichment Analysis
* Slides - gene ontology background and resources.
* Interactive GO exercises and solutions from http://gohandbook.org
* The entire book or just chapters are available online to download
* More code can be found at https://github.com/tanghaibao/goatools/tree/master/notebooks

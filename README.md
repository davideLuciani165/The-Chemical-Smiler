# THE CHEMICAL SMILER
 
**A KNIME Workflow for Context-Dependent Chemical Structure Standardization and Abstraction**
 
## Overview
The **CHEMICAL SMILER** is an open-source, flexible cheminformatics workflow developed within the KNIME Analytics Platform. It bridges the gap between rigid, purely syntactic SMILES canonicalization and the practical needs of predictive toxicology and regulatory meta-analyses.

By separating chemical structures into distinct, customizable semantic layers, the tool allows users to choose exactly which molecular features to preserve (e.g., ionic forms, specific tautomers, components) and which to abstract away based on the specific domain problem.

---

## Prerequisites & Installation

To run the CHEMICAL SMILER, you need to:

### 1. the **KNIME Analytics Platform** (version 4.7 or higher recommended) 
Download and install the core platform from the official website: [KNIME Downloads](https://www.knime.com/downloads).

### 2. KNIME Extensions
Open KNIME and navigate to **File > Install KNIME Extensions...**, then search for and install the following nodes:
* **RDKit KNIME Integration**: Provides robust, standard-compliant molecular canonicalization and structural handling.
* **KNIME Interactive R Statistics Integration**: Required to execute the underlying custom R code embedded within the workflow nodes.

### 3. R
Download and install R (https://www.r-project.org/)

### 4. The R package 'Rserve' 
Open R and install the package Rserve by typing the code: **>install.packages("Rserve")**. 
If you use RStudio, you can click on **Tools > Install packages** and choose the following options:
	- Install from: Repository (CRAN)
	- Packages: Rserve
	- Install dependencies: yes
Once installed, you can close R. To use the CHEMICAL SMILER, you don't need to open R again, only KNIME.
 
### 5. R Environment Setup
Ensure that **R** is installed on your local machine and properly linked to KNIME (**File > Preferences > KNIME > R**). If you need more memory in KNIME to run R, you can increase the Rserve receiving buffer size limit.
 
### 6. Set the available software sites in KNIME
Set the available software sites by clicking on **File > Preferences > Install/Update > Available Software Sites** and check all the boxes; in particular, be sure to check “KNIME Community Extensions”.

### 7. Import the CHEMICAL SMILER
You can import the CHEMICAL SMILER workflow by clicking on **File > Import KNIME workflow**. To open the workflow, double-click on it.

---
## Usage Instructions
1. Importing the Workflow
Download or clone this repository to your local machine.
 
Open KNIME Analytics Platform.
 
Go to **File > Import KNIME Workflow**, select the unpacked root folder, and double-click to start the workflow.
 
2. Configuring the Input Modules and selecting the Abstraction Options
See [The Chemical SMILERS usage](https://github.com/ChemicalSMILER_v1.5.2/TheChemicalSMILERS_HowToUse_mod22July2024.pdf)

3. Running the Workflow and Exporting Results
Right-click the final node and select Execute All.
 
Once execution is complete, right-click the output node to view the newly generated standardized columns or look at the Results folder that contains the output as .csv files.
 
---
## Additional tips
We recommend viewing the heap status and freeing it between runs. To do so, click **File > Preferences > General** and select Show Heap Status.

To free up memory, kill the Rserve process in Task Manager (only if you've completed the run and want to start a new one from the beginning; don't do this if you only want to restart the last "Structure Manipulation" component).


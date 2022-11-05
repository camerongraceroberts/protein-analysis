# Helpful Functions for Protein Analysis

This python program is for calculating protein features from amino acid sequence and extracting information for a given protein from the Uniprot database. The inputs are protein sequence or accession number as a string. 

**You can:**
1. Calculate the number of amino acids, molecular weight in kDa, and the extinction coefficient (reduced or oxidized)
2. Calculate protein concentration in a sample from optical density at 280 nm. 
3. Retrieve a custom protein purification protocol based on protein features such as instability index, isoelectric point, and cysteine content. The protocol is specifically meant for soluble proteins expressed in E. coli.
4. Extract annotation and function data from the Uniprot database. 


## Usage ##
---
The dependencies for this program were managed using [poetry](https://python-poetry.org/)


To install all dependencies: 
```
poetry install 
```


## Using in Python ##
---

**Using the virutal environment:**
```
poetry shell
python -i analysis.py
```

**To determine # of amino acids, MW (kDa), and extinction coefficient from protein sequence**

**Input:** protein sequence as a string 

**Output:** a list [# of Amino acids, Molecular Weight in kDa, extinction coefficient if all Cys are reduced, extinction coefficient if Cys are oxidized]

```
protein_features(sequence)
```


**To calculate protein concentration in sample from optical density**

**Input:** 
1. Optical density of sample at 280 nm
2. If reducing agent in sample input is "y", otherwise "n" 
3. protein sequence as a string 

**Output:** a list [mol/L, mg/mL]

```
calculate_protein_conc(OD280, reducing_agent, sequence)
```

**To retrieve recommended purification protocol for a soluble protein expressed in E. coli**

**Input:** protein sequence as a string 

**Output:** download buffer conditions, materials, and a written protocol (.txt)

```
protein_prep_protocol(sequence, protein_name)
```

**To query the Uniprot database for a protein of interest**

**Input:** protein accession number as a string

**Output:** 
1. Entry Name
2. Gene Name
3. Organism
4. Catalytic activity
5. Protein families
6. Motif search
7. Domains
8. AlphaFold database #

```
download_uniprot_search(accession)
```

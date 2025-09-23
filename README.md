![PolyLattice Banner](https://github.com/andresordorica/Polylattice/blob/main/images/polylattice_banner.png?raw=true)
# PolyLattice
A Python-based toolkit that leverages the functionality of mBuild, part of the MoSDeF simulation suite, to create reproducible crosslinked polymer structures. The toolkit enables high control over the final structure, facilitating visualization and generation of molecular dynamics (MD) simulation inputs that are Transparent, Reproducible, Usable by Others, and Extensible (TRUE; DOI: https://doi.org/10.1080/00268976.2020.1742938).

## Getting Started

To get started, first install the necessary packages to run the tutorials.  
You can create the conda environment from the provided `.yml` file and then activate it:

```bash
conda env create -f conda_environment/original_env.yml -n Polylattice-MoSDeF
conda activate Polylattice-MoSDeF
```

After creating the environment, go into the directory to install the dependencies:

```bash
cd /Polylattice
python -m pip install -e
```


Once the environment and dependencies are ready, begin by exploring the Polymer_Chain_Creation.ipynb notebook
to get familiarized with the MoSDeF–mBuild concepts behind the Polymer class and the
subsequent creation of crosslinked polymer compounds.

This repository contains Examples for:

a) Polymer_Chain_Creation.ipynb: Building polymer compounds from custom repeat units

b) PVA.ipynb:
    PVA crosslinked with Glutaraldehyde, Dimethylolurea and Malic Acid

c) BADGE_MHHPA.ipynb: Crosslinkinging of Bisphenol-A-diglycidyl-ether (BADGE:https://pubchem.ncbi.nlm.nih.gov/compound/Bisphenol-A-diglycidyl-ether)
    and crosslinker Methylhexahydrophtalic-anyhidride (MHHPA: https://pubchem.ncbi.nlm.nih.gov/compound/Methylhexahydrophthalic-anhydride )

d) EPON-862_DETDA.ipynb:
    Crosslinking of a polyamide EPON 862 with the curing Agent DETDA (https://pubs.acs.org/doi/full/10.1021/ma2005519)

e) TETA-TMC*.ipynb:
    Crosslinking of a polyamide triethylentetramine (TETA) with trimesoyl chloride (TMC) (https://pubs.acs.org/doi/10.1021/acsomega.3c10108)
    I)  TETA-TMC_via_API.ipynb: Using PubChem's API for loading the initial monomers
    II) TETA-TMC_via_MOL2.ipynb: Using mol2 files 

The examples shown in this repository are ilustrative as the functionality of the code can be adapted to many other chemistries using the
methodology and outlined strategies shown in the current repository

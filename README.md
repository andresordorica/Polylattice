![PolyLattice Banner](https://github.com/andresordorica/Polylattice/blob/main/images/polylattice_banner.png?raw=true)
# PolyLattice
A Python-based toolkit that leverages the functionality of mBuild, part of the MoSDeF simulation suite, to create reproducible crosslinked polymer structures. The toolkit enables high control over the final structure, facilitating visualization and generation of molecular dynamics (MD) simulation inputs with a strong emphasis on reproducibility.

## Getting Started

To get started, first install the necessary packages to run the tutorials.  
You can create the conda environment from the provided `.yml` file and then activate it:

```bash
conda env create -f conda_environment/original_env.yml -n Polylattice-MoSDeF
conda activate Polylattice-MoSDeF
```


Once the environment is ready, begin by exploring the Creating_polymer.ipynb notebook
to get familiarized with the MoSDeF–mBuild concepts behind the Polymer class and the
subsequent creation of crosslinked polymer compounds.

This repository contains examples for:

a) PVA.ipynb:
    PVA crosslinked with Glutaraldehyde and Dimethylolurea

b) BADGE_MHHPA.ipynb:

c) EPON-862_DETDA.ipynb:
    Crosslinking of a polyamide EPON 862 with the curing Agent DETDA (https://pubs.acs.org/doi/full/10.1021/ma2005519)

c) EPON-862_DETDA.ipynb:
    Crosslinking of a polyamide EPON 862 with the curing Agent DETDA (https://pubs.acs.org/doi/full/10.1021/ma2005519)

The examples shown in this repository are ilustrative as the functionality of the code can be adapted to many other chemistries using the
methodology and outlined strategies shown in the current repository

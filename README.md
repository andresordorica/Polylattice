![PolyLattice Banner](https://github.com/andresordorica/Polylattice/blob/main/images/polylattice_banner.png?raw=true)
# PolyLattice
A Python-based toolkit that leverages the functionality of mBuild, part of the MoSDeF simulation suite, to create reproducible crosslinked polymer structures. The toolkit enables high control over the final structure, facilitating visualization and generation of molecular dynamics (MD) simulation inputs with a strong emphasis on reproducibility.

## Getting Started

To get started, first install the necessary packages to run the tutorials.  
You can create the conda environment from the provided `.yml` file and then activate it:

```bash
conda env create -f conda_environment/original_env.yml -n Polylattice-MoSDeF
conda activate Polylattice-MoSDeF



Once the environment is ready, begin by exploring the Creating_polymer.ipynb notebook to get familiarized with the MoSDeF–mBuild concepts behind the Polymer class and the subsequent creation of polymer compounds.
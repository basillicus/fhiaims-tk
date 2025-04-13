# Traincraft Toolkit

Set of independent python scripts to facilitate simple but frequent tasks you may need to perform while working with electronic structure calculation files.

These little scripts perform operations that may be needed when manipulating molecular geometries, inverconverting formats or handling keywords/units between different codes. 

This is a complement to the code Traincraft, for the generation of datasets for traininf Machine Learning Interatomic Potentials.

It uses ASE as external dependency, which depends on numpy, and some scripts need scikit-learn

## Configure

1. Clone the repository.
2. Copy the config_template.py file to config.py file and edit it to add your own paths. (not all scripts make use of it)
3. You need a python interpreter configured in your system and ASE installed. I like to have a dedicated conda environment:

       conda create -c conda-forge -n compchem python=3 ase scikit-learn
       conda activate compchem

You may need also tqdm if you use some of the scripts that use it. Some files are very large and takes a while for them to read and parse them, that is why it is better to see a progress bar. If you need tqdm:

    (compchem)$ conda install tqdm


### How to use it

I create an alias that adds to $PATH the path where the scripts are stored.

In ~/.bash_aliases:

    alias load_tctk_tools="export PATH=$PATH:$HOME/git/traincraft-tk/"

When you need the scripts, activate your environment, and load the alias:

    conda activate compchem
    load_tctk_tools

All scripts are executable and start with "tctk_", so now you can execute each individual script. You can pass options to the scripts. For more information type for example:

    tctk_<name-of-the-script>.py -h

# Acknowledgment
This project has received funding from the European Union’s Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie grant agreement No 847635, under the project UNA Europa, an alliance of universities FOR the emergence of talent and the development of research CAREERs (UNA4CAREER)

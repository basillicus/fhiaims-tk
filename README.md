# FHI-AIMS Toolkit

This is a very simple set of python scripts to facilitate simple but frequent tasks you may need to perform while working with FHI-AIMS
I've just started working with FHI-AIMS and probably this tools already exist somewhere else, or you do it in a different and probably more efficient way.

At the moment there are no many tools, I guess I will keep adding them as long as I keep needing them and can not find them.

At some point it will use ASE, but at the moment does not require any external dependence.


## Configure

1. Clone the repository
2. Edit the config.py file to add your own paths
3. You need a python interpreter configured in your system. I like to have a dedicated conda environment:

    conda create -c conda-forge -n compchem python=3
    conda activate compchem

### My own set up  

I do the following because I believe is simple, and does not mess with your path until you need it

I create an alias that adds to $PATH the path where the scripts are stored

In ~/.bash_aliases:

    alias load_fhi_tools="export PATH=$PATH:$HOME/git/fhiaims-tk/"

All scripts are executable and start with "fhi_", so now you can execute each individual script.To some scripts you can pass options to them. For example:

    fhi_add_basisSet.py -h


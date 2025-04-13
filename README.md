# Traincraft Toolkit

Set of independent python scripts to facilitate simple but frequent tasks you may need to perform while working with electronic structure calculation files
I've just started working with FHI-AIMS and probably this tools already exist somewhere else, or you do it in a different and probably more efficient way.

At the moment there are no many tools, I guess I will keep adding them as long as I keep needing them and can not find them.

It uses ASE as external dependency, which depends on numpy, and some scripts need scikit-learn


## Configure

1. Clone the repository.
2. Copy the config_template.py file to config.py file and edit it to add your own paths. (not all scripts make use of it)
3. You need a python interpreter configured in your system and ASE installed. I like to have a dedicated conda environment:

       conda create -c conda-forge -n compchem python=3 ase scikit-learn
       conda activate compchem

You may need also tqdm if you use some of the scripts that use it. Some files are very large and takes a while for them to read and parse them, that is why it is better to see a progress bar. If you need tqdm:

    (compchem)$ conda install tqdm


### My own set up  

I do the following because I believe it is simple, and does not mess with your path until you need it.

I create an alias that adds to $PATH the path where the scripts are stored.

In ~/.bash_aliases:

    alias load_tctk_tools="export PATH=$PATH:$HOME/git/traincraft-tk/"

When you need the scripts, activate your environment, and load the alias:

    conda activate compchem
    load_tctk_tools

All scripts are executable and start with "tctk_", so now you can execute each individual script. You can pass options to the scripts. For more information type for example:

    tctk_<name-of-the-script>.py -h


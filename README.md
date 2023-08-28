# The locmodslib a small library for Local Mode Approach #

Collection of a small library and few scripts to perform a local mode calculation.

## How to ##

Requires a Gaussian log file where read geometry and few data (es. gaus.log) and a json file with Gaussian16 parameters (Functional, basis sets...).

### Install ###
From the root folder activate or create a Python virtual environments and run:

`pip install -e .`

#### Requirements ####

 - `numpy`, `scipy`, [`estampes`](https://github.com/jbloino/estampes)

### Steps ###

- Generate the structure for the selected XH bonds

`lmmkinp gaus.log opts.json "C" -w test --prefix gaus`

- Run all the gaussian input files, in utilities an example of bash script is provided which crops the fchk in order to reduce the disk usage.
On a node run:

`bash runlmodes.sh gaus 2`

- Collect the data, all the gjf files are expected in the same folder with the types of name written by `lmodesmkinp.py`

`lmprcdt gauss.log "C" "test" "gaus"`

### Build the Docs ###

Requires sphinx, sphinx-argparse, sphinx-rtd-theme

`cd docs; make html`





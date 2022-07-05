# locmodslib a small library for Local Mode Approach #

Collection of a small library and few scripts to perform a local mode calculation.

## How to ##

Requires a Gaussian log file where read geometry and few data (es. gaus.log) and a json file with Gaussian16 parameters (Functional, basis sets...).

### Steps ###

- Generate the structure for the selected XH bonds

`python lmodesmkinp.py gaus.log opts.json "C" -w test --prefix gaus`

- Run all the gaussian input files, in utilities an example of bash script is provided which crops the fchk in order to reduce the disk usage.
On a node run:

`bash runlmodes.sh gaus 2 `

- Collect the data, all the gjf files are expected in the same folder with the types of name written by `lmodesmkinp.py`

`python processjobs.py gauss.log "C" "test" "gaus"`





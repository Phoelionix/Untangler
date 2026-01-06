# Untangler

Altloc optimizer to [untangle confused conformers](https://bl831.als.lbl.gov/~jamesh/challenge/twoconf/).

## Setup 
(Optional) make a virtual Python environment

`python3.9 -m venv untangler_venv`

`source untangler_venv/bin/activate`

Install required packages

`pip install -r requirements.txt` 

Tested with python 3.9. Requires `phenix`.

## Run demo

`python3.9 untangle.py data/longrangetraps_TW.pdb data/refme.mtz`


<!-- <CENTER><P>
<HR><A href="untangling.gif"><img src=untangling.gif width=960 height=720></A><p>
<HR></P></CENTER> -->


## What is this?


This program refines an ensemble of model protein conformations to fit the 
electron density represented by X-ray data (using `phenix.refine`), while 
simultaneously freeing the models from local minima traps that occur when a 
the conformers of a single model conformation are "fitting" the electron density 
of multiple true conformations. The goal is to reallot or "reconnect" conformers 
in a way that not only improves the fit, but also reduces the potential energy 
of the conformations.

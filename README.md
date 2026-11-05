# Untangler

Connectivity optimization for reducing strain in ensemble models.

This readme is out of date!

## Requirements
`phenix`, `python3.12` or later, `tcsh`, `gnuplot`, 

## Setup 
(Optional) make a virtual Python environment

`python3 -m venv untangler_venv`

`source untangler_venv/bin/activate`

Install required packages

`pip install -r requirements.txt` 

Tested with python 3.12.

## Untangle Challenge example 

Application to synthetic dataset of the [Untangle Challenge](https://doi.org/10.64898/2026.02.21.706873).

IMPORTANT NOTE: On this branch, the program will NOT arrive at the  solution on its default settings, which have been chosen so to improve the connectivity of experimental models with more than two conformations. Notably, untwist moves are disabled, a single suboptimal solution is identified rather than a batch of highest-scoring solutions, and only the weighted energy of a single conformation of each model is tracked (pending the implementation of a parallel computation of conformation geometry scores). 
For a demonstration, please switch to the branch: https://github.com/Phoelionix/Untangler/tree/2_conformer_challenge_solution


`python3 untangle.py data/synthetic/longrangetraps.pdb data/synthetic/refme.mtz`  

To automatically generate information on the differences from the best model run

`python3 untangle.py data/synthetic/longrangetraps.pdb data/synthetic/refme.mtz data/synthetic/best.pdb`  

This will generate a file like `output/linear_optimizer_logs/cheat_longrangetraps_[...]/ChangedBonds-1.txt`. 


<!-- <CENTER><P>
<HR><A href="untangling.gif"><img src=untangling.gif width=960 height=720></A><p>
<HR></P></CENTER> -->


## What is this?

This program refines an ensemble model of protein conformations to fit the 
electron density represented by X-ray data (using `phenix.refine`), while 
simultaneously freeing the models from local minima traps that occur when 
the conformers of a single model conformation are "fitting" the electron density 
of multiple true conformations. The goal is to reallot conformer labels to "reconnect"
them in conformations with better geometries and which better fit the X-ray data.

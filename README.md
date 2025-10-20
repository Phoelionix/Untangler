# Untangler

Code for treating [confused conformers](https://bl831.als.lbl.gov/~jamesh/challenge/twoconf/).

## Setup 
(Optional) make a virtual Python environment

`python3.9 -m venv untangler_venv`

`source untangler_venv/bin/activate`

Install required packages

`pip install -r requirements.txt` 

Tested with python 3.9. Requires `phenix`.

## Run demo

`python3.9 untangle.py data/longrangetraps_TW.pdb data/refme.mtz`


<CENTER><P>
<HR><A href="untangling.gif"><img src=untangling.gif width=960 height=720></A><p>
<HR></P></CENTER>


## What is this?


This program refines an ensemble of model protein conformations to fit the 
electron density represented by X-ray data (using `phenix.refine`), while 
simultaneously freeing the models from local minima traps that occur when a 
the conformers of a single model conformation are "fitting" the electron density 
of multiple true conformations. The goal is to reallot or "reconnect" conformers 
in a way that not only improves the fit, but also reduces the potential energy 
of the conformations.


<!-- ## Background 
### Ensemble modelling
Proteins are big and wiggly. In a crystal, they can take up multiple geometric conformations.
This means that the data we collect does not represent the electron density of a single structure,
but rather, it is the overlaid densities of many different geometric conformations.

This poses an issue: the data does not tell us the contribution by each conformation 
to the total electron density at a given point in space. 
However, we can get a very good idea of how bad our model is by considering the geometric
'badness' of the structure. If a conformation 'wants' to relax to a different conformation,
then we can be pretty certain it wasn't a conformation within our protein crystal.

  

### The problem

Necessitating that the electron density is fit with only relaxed geometric structures is a 
powerful constraint. Indeed, this is also the strategy of single-conformation model refinement, 
where atoms are moved to improve the fit while maintaining a geometrically plausible structure. 
However, this is done through many small steps. 

Unfortunately, refining ensembles in this way even against noiseless data 
results in models that are trapped in local minima, where improving one 
measure (geometry or X-ray) will make the other worse*. 
Specifically, these local minima correspond to cases where a model conformation contains, 
and thus bonds, two atoms that are fitting the electron density of two different structures. 
In these cases, we need to make big changes. 

*It's worth taking a moment to appreciate that this in itself indicates an incorrect model.
We expect both measures to be at or near local minimums. -->




<!-- Want to prevent models from fitting the electron density of multiple conformations.
This code tries to find the 'correct' way to connect atoms in ensemble models. -->
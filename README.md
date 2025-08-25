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

NB: Currently the code assumes the input .pdb and .mtz files are in `data/`


<CENTER><P>
<HR><A href="untangling.gif"><img src=untangling.gif width=960 height=720></A><p>
<HR></P></CENTER>

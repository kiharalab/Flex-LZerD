# Flex-LZerD

Flex-LZerD models the bound conformation of a protein complex, where one protein binds to another macromolecule after undergoing an extreme conformational change. It can also be used to morph one protein conformation into another, even if the structure of the endpoint is incomplete.

Copyright (C) 2022 Charles Christoffer, Daisuke Kihara, and Purdue University.

License: GPL v3

Contact: Daisuke Kihara (dkihara@purdue.edu)

## Citation:
```
@article{cwc2022flex,
  title={},
  author={Charles Christoffer, Daisuke Kihara},
  journal={Submitted for publication}
}
```

Installation
============

To run the test, all Python dependencies are required.

System requirements
-------------------
All Flex-LZerD tools were developed and tested on the CentOS and Ubuntu families of Linux operating systems. Flex-LZerD does not actively support Windows, MacOS, or other environments. The minimum suggested hardware is an AMD Rome or Intel Skylake-SP or newer CPU and at least 8 GB RAM.

Python dependencies
-------------------
- python2.7
- Biopython
- numpy
- scipy
- prody

In an anaconda environment, you can run `conda install -c insilichem prody python=2.7 numpy scipy Biopython` to install these Python dependencies.

Binary dependencies
-------------------
- LZerD ([https://kiharalab.org/proteindocking/lzerd.php](https://kiharalab.org/proteindocking/lzerd.php))
- GOAP ([http://cssb.biology.gatech.edu/GOAP/index.html](https://sites.gatech.edu/cssb/goap/))
- ITScorePro ([http://zoulab.dalton.missouri.edu/resources.html](http://zoulab.dalton.missouri.edu/resources_itscorepro.html))
- Phenix ([https://phenix-online.org](https://phenix-online.org)) version 1.17 or newer.

Getting started
===============
Before running any instance of flexible fitting, the Phenix suite must be loaded into your environment. Instructions for activating your installation of Phenix are included in the [Phenix documentation](https://phenix-online.org/documentation/install-setup-run.html#setting-up-the-command-line-environment). If your environment is managed via Lmod or other module system, the proper way to load Phenix may differ. In such a case, you should ask your system administrator about your system-specific way to load Phenix. One common Lmod command to load Phenix is `module load phenix`.

To run flexible fitting on a structure, you should first determine your unbound ligand structure and domains suitable for rigid-body docking. Extracted domains should be reasonably contiguous, that is, they should not cut the protein sequence many times. Flex-LZerD will not reject such inputs and in principle can function using them, but strange and extreme extractions such as cutting a globular region in half without regard for the topology can interfere with docking.

Test protein
------------
Flexible fitting can be run for the provided example protein and domains found in the `example/` directory using the command:

```run_flexible_fitting.py docked_domains.pdb input_ligand_from_3nd2.pdb x x x gg1.0myout receptor_from_3ea5.pdb```

You should allow several hours for flexible fitting to complete. At each iteration, an output structure is written to `lastiter.pdb`, which can be viewed at any point during the fitting. At the end, your output should be similar to the example output provided in `example/flex_output_ligand_model.pdb`.

Running a new protein
---------------------
Once you have extracted your domains as above, you can then dock them using a protein docking method, for example LZerD. LZerD is available via an interactive web server at [https://lzerd.kiharalab.org](https://lzerd.kiharalab.org) and via downloadable executable at [https://kiharalab.org/proteindocking/lzerd.php](https://kiharalab.org/proteindocking/lzerd.php). Once you have obtained the domains you want to fit to, simply combine them into one PDB file, e.g. via [PyMOL](https://pymol.org) or the molecular editor of your choice, and follow the command as above.

Output interpretation
---------------------
Flex-LZerD will generate a `.tar` file containing all frames of the flexible fitting to your domains, as well as a `lastiter.pdb` file containing the most recent frame. By unzipping the tarball and loading the contained PDB frame files into PyMOL, you can see the progression of your ligand as it is fit to your domains. An example final frame is given in `example/flex_output_ligand_model.pdb`.

Scoring
-------
The residue occupancies as described in the paper can be calculated via
```
count_residue_decoy_contacts receptor.pdb ligand.pdb ${RADIUS} < lzerd-rigid-transformations.txt > occupancies_${RADIUS}.txt'
```
where RADIUS is the occupancy distance cutoff, default 5.0 Å. The binding site consensus scores can then be calculated via
```
calc_binding_consensus_score receptor.pdb ligand.pdb ${RADIUS} occupancies_${RADIUS}.txt < lzerd-rigid-transformations.txt > binding_site_consensus_${RADIUS}.txt
```

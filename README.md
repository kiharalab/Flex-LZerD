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

Flexible fitting can be run for the provided example protein and domains using the command:

```run_flexible_fitting.py docked_domains.pdb input_ligand_from_3nd2.pdb x x x gg1.0myout receptor_from_3ea5.pdb```

You should allow several hours for flexible fitting to complete. At each iteration, an output structure is written to `lastiter.pdb`, which can be viewed at any point during the fitting.

Running a new protein
---------------------

Scoring
---------------------

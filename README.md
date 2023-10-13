# Molecule_alignment
Automated molecule alignment script from 2D structure sdf files. The 3D coordinates of XYZ files will be generated.

# Script Usage

Prepare a `<mols>.sdf` file containing all molecules to be aligned and `<ref>.sdf` as an alignment reference.

Run the command : `python3 align.py --ref <ref>.sdf --mols <mols>.sdf --type mmff`

# Reference

[https://iwatobipen.wordpress.com/2018/09/26/3d-alignment-function-of-rdkit-rdkit/](https://iwatobipen.wordpress.com/2018/09/26/3d-alignment-function-of-rdkit-rdkit/)

# Perovskite_Analysis
1) Extract structural features of 2D, 3D, and 2D/3D double perovskite systems.
2) Build and analyze bulk perovskites, double perovskites, 2D perovskites, and 2D double perovskites.
3) Generate elctronic structure input files for xTB/GPAW

## Requirements
- ASE https://wiki.fysik.dtu.dk/ase/
- pymatgen https://pymatgen.org/
- pytest to run ```pytest .``` in the tests directory.

## Installation
The package can be installed via 
```python
pip install pyrovksite
```

or you can clone the repository locally with:

```
git clone https://github.com/r2stanton/pyrovskite.git .
```

With this method you'll download the tests and example jupyter notebooks. Tests can be run with pytest by navigating to the tests/ folder and running pytest.


## Note on xTB input file generation
xTB input files are for use with CP2K. Additionally, openbabel is required to
convert files to a CP2K readable format. You can install both of these things
with the command.
 - cp2k https://www.cp2k.org
 - openbabel http://openbabel.org/wiki/Main_Page
 
```python
conda install -c conda-forge cp2k openbabel
```

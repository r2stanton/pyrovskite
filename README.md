# Pyrovskite
1) Extract structural features of 2D, 3D, and 2D/3D double perovskite systems.
2) Build and analyze bulk perovskites, double perovskites, 2D perovskites, and 2D double perovskites.
3) Generate elcetronic structure input files for xTB/GPAW.

#### 🟩 If you use this code for your research, please cite: https://doi.org/10.1063/5.0159407 🟩

## Requirements
- ASE https://wiki.fysik.dtu.dk/ase/
- pymatgen https://pymatgen.org/
- pytest to run ```pytest .``` in the tests directory.

## Installation
The package can be installed via 
```python
pip install pyrovksite
```

or you can clone the repository locally and obtain an editable install with:

```
git clone https://github.com/r2stanton/pyrovskite.git .
pip install -e .
```

With this method you'll download the tests and example jupyter notebooks. Tests can be run with pytest by navigating to the `tests/` folder and running ```pytest .```

A small collection of organic A'-site and A-site are available in the `res/` folder, see `spacers.csv` for charge info.

## Usage
A brief overview of code usage is discussed here, however a much more detailed explanation is contained in the Jupyter notebooks found in `/examples`.

### Example usage, from structure creation to analysis.

Example of building a bulk MAPbI3 perovskite:

```python
from pyrovskite.builder import make_bulk
import ase.io
A = ase.io.read("path/to/methylammonium.xyz")
MAPbI3 = make_bulk(A, "Pb", "I", 3.1)
```
Similar such methods exist for 2D phases, and double perovskites. These are returned as an Atoms object, allowing for all the I/O methods supported by ASE.

One can write an input file, for example to do a variable cell xTB optimization with:
```python
MAPbI3.write_xTB()
```

The optimized geometry could then be obtained and loaded back into a Perovksite object for further analysis:
```python
from pyrovskite.perovskite import Perovskite
MAPbI3 = Perovskite("path/to/optimized_MAPbI3.cif")
```

Structural properties such as octahedral distortions can be computed with ```MAPbI3.compute_delta()``` ```MAPbI3.compute_sigma()``` ```MAPbI3.compute_lambda()``` and plots of relevant distributional quantities such as bond angle, bond length, and pRDFs can be computed with ```MAPbI3.plot_angles()``` ```MAPbI3.plot_distances()``` ```MAPbI3.plot_rdf()```

Sensible defaults are selected for these functions, but please refer to the in-built documentation with ```help(function)``` or to the example notebooks to ensure the default behavior is what you want.


### Note on xTB input file generation
xTB input files are for use with CP2K. Additionally, openbabel is required to
convert files to a CP2K readable format. You can install both of these things
with the command.
 
```python
conda install -c conda-forge cp2k openbabel
```
 - cp2k https://www.cp2k.org
 - openbabel http://openbabel.org/wiki/Main_Page

Additionally, if you plan on using the xTB optimizations frequently, it may be worth setting 
```zsh
export XTB_D3_PATH=/path/to/dftd3.dat
export XTB_PARAMS_PATH=/path/to/xTB_parameters
```
to your shell profile, as the xTB calculation requires these parameters. If they are declared as environment variables, then they will be automatically parsed for you, otherwise you need to make a custom input file template, or pass the paths to these parameter files as keyword arguments to the relevant functions.

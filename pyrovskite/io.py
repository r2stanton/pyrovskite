from pymatgen.io.ase import AseAtomsAdaptor
import ase.io, os

def xTB_input(calc_type = "vcopt", prefix = "xtb", cif_name = None,
              atoms = None, xtb_path = None, d3_path = None):

    """
    Generate input file for xTB calculation (as implemented in CP2K).

    You must either set XTB_D3_PATH, and XTB_PARAMS_PATH as environment
    variables, or pass them into the relevant keywords (see below).

    Recommended to add e.g.
    export XTB_D3_PATH=/path/to/dftd3.dat
    export XTB_PARAMS_PATH=/path/to/xTB_parameters
    to your shell profile. If you do this, python will parse your 
    environment variables for these file paths, and passing to the
    funtion is not needed.

    This function has two modes:
    ---------------------------MODE 1--------------------------------------
    calc_type = "custom" | Where you bring your own template file, and set this
    in the template kwarg. You can use the replace dictionary to store patterns
    for replacement with information about your system.

    E.g. replace = {"@halogen_corr":".TRUE."} will search for @halogen_corr
    within your template file, and replace it with .TRUE. when generating
    the input file. 
    ---------------------------MODE 1--------------------------------------

    ---------------------------MODE 2--------------------------------------
    calc_type = "vcopt", "opt", "md", "scf" | Where you use one of the
    default files contained explicitly, line for line in this function. 

    The template used is determined by calc_type. The only things that will
    be determined on the fly are UKS (if nelect is odd) and the cif_file used
    to provide CP2K with cell parameters and atomic coordinates.

    Feel free to take these templates, modify them manually, and swap to 
    calc_type = "custom" mode for more control over the templating. In fact
    this is encouraged, as the default templates are OVERLY general, and
    certainly not optimized for any particular system.
    ---------------------------MODE 2--------------------------------------

    Parameters:
    ----------
    calc_type (str)
        Type of calculation. Default is "vcopt". This
        is currently also the only option aside from "custom". See above.

    prefix (str)
        Prefix for output files. Default is "xtb". If a cif_name is not
        provided, but an Atoms object is, a cif of the form:
        f'{prefix}_{calc_type}.cif' will be generated and used to feed to
        CP2K.

    atoms (ase.Atoms)
        Atoms object to be used to generate cif file. If this is not
        provided, a cif_name must be provided. Don't provide cif_name
        if atoms are being passed.

    cif_name (str)
        Name of cif file to be used to feed to CP2K. If this is not
        provided, an atoms must be provided. Don't provide
        atoms if cif_name is being passed.
        
    Returns:
    -------
        None, generates input file:
        f'{prefix}_{calc_type}.inp'
    """

    assert cif_name is not None or atoms is not None, \
        "Must provide either cif_name or atoms"

    assert cif_name is None or atoms is None, \
        "Must provide either cif_name or atoms, not both"


    # If no cif, needs to be written.
    if atoms is not None:
        cif_name = f'{prefix}_{calc_type}.cif'
        atoms.write(cif_name)

        # Let's try to use pymatgen to write the cif file.
        # Doesn't work
        # AseAtomsAdaptor().get_structure(atoms).to(filename=cif_name,
                                                  # fmt="cif")
        os.system(f"obabel {cif_name} -O {cif_name}")

        # Determins if UKS is needed.
        if _count_electrons(atoms) % 2 == 1:
            UKS = ".TRUE."
        else:
            UKS = ".FALSE."
    else:
        if _count_electrons(cif_name) % 2 == 1:
            UKS = ".TRUE."
        else:
            UKS = ".FALSE."


    # Checks kwargs
    if xtb_path is None:
        xtb_err_string = "Must provide path to xTB parameters file either" \
                         "as an environment variable, or to the function.\n" \
                         "See help(xTB_input) for more information."


        # Checks environment variables, raises error if unset.
        xtb_path = os.getenv("XTB_PARAMS_PATH")
        if xtb_path is None:
            raise ValueError(xtb_err_string)
        else:
            if type(xtb_path) == str:
                if len(xtb_path) == 0:
                    raise ValueError(xtb_err_string)


    # Checks kwargs
    if d3_path is None:
        d3_err_string = "Must provide path to D3 parameters file either" \
                         "as an environment variable, or to the function.\n" \
                         "See help(xTB_input) for more information."

        # Checks environment variables, raises error if unset.
        d3_path = os.getenv("XTB_D3_PATH")
        if d3_path is None:
            raise ValueError(d3_err_string)
        else:
            if type(d3_path) == str:
                if len(d3_path) == 0:
                    raise ValueError(d3_err_string)



    if calc_type == "vcopt":
        with open(f"{prefix}_vcopt.inp", "w") as fil:
            fil.write(f"""\
&FORCE_EVAL
  &DFT
	UKS {UKS}
    &QS
      METHOD xTB
      &xTB
        DO_EWALD  T
        CHECK_ATOMIC_CHARGES  F
        COULOMB_INTERACTION T
        &PARAMETER
          DISPERSION_PARAMETER_FILE {d3_path}
          PARAM_FILE_NAME {xtb_path}
        &END PARAMETER
        USE_HALOGEN_CORRECTION .TRUE.
      &END

      &DISTRIBUTION
        BASIC_OPTIMIZATION .FALSE.
        BASIC_SPATIAL_OPTIMIZATION .TRUE.
      &END
      
    &END QS
    &POISSON
      POISSON_SOLVER PERIODIC
      PERIODIC XYZ
    &END

    &SCF
      SCF_GUESS ATOMIC
      EPS_SCF 1.e-8
      &OT
         PRECONDITIONER FULL_SINGLE_INVERSE
         MINIMIZER DIIS
		 ENERGY_GAP .03 ! A conservative ~.8 eV for most perovskites.
      &END
      MAX_SCF  500
    &END SCF

  &END DFT

  STRESS_TENSOR ANALYTICAL

  &SUBSYS
    &CELL
        CELL_FILE_FORMAT CIF
        CELL_FILE_NAME {cif_name}
    &END CELL
    &TOPOLOGY
        COORD_FILE_NAME {cif_name}
        COORD_FILE_FORMAT CIF

      &GENERATE 
        !REORDER T 
      &END GENERATE

    &END TOPOLOGY
  &END SUBSYS

&END FORCE_EVAL

&GLOBAL
  PROJECT_NAME {prefix}_{calc_type}
  RUN_TYPE CELL_OPT
  PRINT_LEVEL MEDIUM
&END GLOBAL

&MOTION
  &CELL_OPT
    MAX_ITER 500
    KEEP_SYMMETRY .TRUE.
	KEEP_ANGLES .TRUE.
  &END
  &GEO_OPT
      OPTIMIZER CG 
      MAX_ITER   5000
      MAX_FORCE  9.7225D-4 
      TYPE MINIMIZATION
  &END
  &PRINT
    &FORCES ON
    &END FORCES
  &END PRINT
&END

                    """)

def GPAW_input(calc_type = "fcopt", prefix = "gpaw", cif_name = None,
              atoms = None, opt_algo = "FIRE", fmax = .02):

    """
    Generates basic GPAW input file. Feel free to use as a template for the
    specific calculations you want to run. There are a lot of things that 
    can speed up this type of calculation which are a bit too specific to
    include in a general template like the one used below.

    Parameters
    ----------
    calc_type : str
        Type of calculation to run. Options are "fcopt" currently.
    prefix : str
        Prefix for all output files.
    cif_name : str
        Name of cif file to use. If None, atoms must be provided.
    atoms : ase.Atoms
        Atoms object to use. If None, cif_name must be provided.
    opt_algo : str
        Optimization algorithm to use. Options are "FIRE" and "BFGS".
    fmax : float
        Force convergence criteria for optimization (in eV/Angstrom).

    Returns
    -------
        None, generates input file:
        f'{prefix}_{calc_type}.py'
    """

    assert cif_name is not None or atoms is not None, \
        "Must provide either cif_name or atoms"

    assert cif_name is None or atoms is None, \
        "Must provide either cif_name or atoms, not both"

    if atoms is not None:
        cif_name = f'{prefix}_{calc_type}.cif'
        atoms.write(cif_name)


    kpt_str = "kpts={'density': 6.0, 'gamma': True},"
    par_str = "parallel={'sl_auto': True},"

    if calc_type == "fcopt":
        with open(f"{prefix}_fcopt.py", "w") as fil:
            fil.write(f"""\

from ase.optimize import BFGS, FIRE
from gpaw import GPAW, PW
import ase.io, sys

curr_structure = ase.io.read('{cif_name}')

calc = GPAW(xc='PBE',
            mode=PW(450, dedecut='estimate'),
            {kpt_str}
            {par_str}
            txt= '{prefix}' + '_fixopt.txt')

curr_structure.calc = calc
fix_relax = BFGS(curr_structure)
fix_relax.run(fmax={fmax})

curr_structure.write('{prefix}' + '_OPT.cif')
                      """)




def _count_electrons(cif):
    """Count the number of electrons in a CIF file or Atoms object.

    Parameters
    ----------

    cif : str or Atoms object
        The CIF file or Atoms object to count electrons of.

    Returns
    -------
    Z : int
        The number of electrons in the CIF file or Atoms object.
    """
    
    if type(cif) == str:
        atoms = ase.io.read(cif)
    else:
        atoms = cif.copy()

    structure = AseAtomsAdaptor.get_structure(atoms)
    chemical_symbols = structure.species
    Z = 0
    for elem in chemical_symbols:
        Z += elem.Z

    return Z


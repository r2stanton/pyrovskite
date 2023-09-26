from evgraf import find_inversion_symmetry
import numpy as np




class Featurizer():
    """
        Featurizer class for applying descriptors to perovkskite data. The
        intent of this class is to be used for featurizing datasets in the form
        of a pandas DataFrame, however the same functionality can be achieved
        by using the featurize method on a single system. 

        The featurizer has two modes.
        -----------------------------------------------------------------------
        'single':
        Featurize a single perovskite system. If this is the case, the user
        must pass the perovskite system as an ASE Atoms object to each method

        For example:
        >>> import ase.io
        >>> my_perovskite = ase.io.read("my_perovskite.cif")
        >>> featurizer = Featurizer()
        >>> featurizer.my_new_feature(perovskite = my_perovskite)

        The empty constructor sets the mode to 'single' automatically.

        -----------------------------------------------------------------------

        -----------------------------------------------------------------------
        'df':
        Featurize a full dataframe of perovskites with the same feature. In this
        case, the user must pass a pandas DataFrame to the constructor, and
        specify the column name which contains the ASE Atoms objects.
        
        For example:
        >>> featurizer = Featurizer(df = my_dataframe, ase_col = 'ASE atoms')
        >>> df = featurizer.my_new_feature()

        This would return a df with the existing data and new feature. Passing
        df and ase_col internally sets the mode to 'df'.

        -----------------------------------------------------------------------

        Note: This approach is intended to be used for features which are not
        'important enough' to be stored within the pyrovskite.perovskite.Perovskite
        objects themselves.

        Also as a temporary note, the 'single' mode is probably not well tested/
        potentially not even implemented for some methods yet.

    """
    def __init__(self, df = None, ase_col = None):
        self.df = df
        self.ase_col = ase_col

        if df is not None and ase_col is not None:
            self.mode = 'df'
        else:
            self.mode = 'single'



    def org_inorg_ratio(self, perovksite = None):

        self._input_consistency_check(perovksite)
        org_syms = ["C", "H", "N", "O"]

        if self.mode == 'single':
            syms = perovskite.get_chemical_symbols()
            org_syms = [sym for sym in syms if sym in org_syms]
            return len(org_syms) / len(syms)

        if self.mode == 'df':
            self.df['org_inorg_ratio'] = self.df[self.ase_col].apply(lambda x: len([sym for sym in x.get_chemical_symbols() if sym in org_syms]) / len(x.get_chemical_symbols()))
            return self.df

    def org_inorg_weight_ratio(self, perovskite = None):
        self._input_consistency_check(perovskite)

        org_syms = ["C", "H", "N", "O"]

        if self.mode == 'single':
            raise NotImplementedError("This method is not yet implemented for single mode.")

        elif self.mode == 'df':

            # Define inner function for application to dataframe rows.
            def get_mass_ratio(atoms):
                syms = atoms.get_chemical_symbols()
                org_indices = [at.index for at in atoms if at.symbol in org_syms]
                org_mass = np.sum(atoms.get_masses()[org_indices])
                total_mass = np.sum(atoms.get_masses())
                return org_mass / total_mass

            self.df['org_inorg_weight_ratio'] = self.df[self.ase_col].apply(get_mass_ratio)

            return self.df

    def measure_inversion_symmetry(self, perovskite = None, progress_apply = True):
        try:
            import tqdm
        except:
            progress_apply = False
            print("TQDM not installed. Progress bars will not be shown.")

        self._input_consistency_check(perovskite)

        def rmsd_from_inversion_symmetry(ats):
            rmsd = find_inversion_symmetry(ats).rmsd
            if rmsd < 1e-5:
                return 0
            else:
                return rmsd
        if progress_apply:
            self.df['rmsd_from_inversion'] = self.df[self.ase_col].progress_apply(rmsd_from_inversion_symmetry)
        else:
            self.df['rmsd_from_inversion'] = self.df[self.ase_col].apply(rmsd_from_inversion_symmetry)

    
    def _input_consistency_check(self, perovskite):
        if self.mode == 'single' and perovskite is None:
            raise ValueError("Must pass a perovskite to each method when in 'single' mode.")
        elif self.mode == 'df' and (self.df is None or self.ase_col is None):
            raise ValueError("Must pass a dataframe and ase_col when in 'df' mode.")





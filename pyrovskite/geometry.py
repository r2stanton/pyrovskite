from ase.neighborlist import natural_cutoffs
from ase.geometry.analysis import Analysis
from ase.visualize import view
from ase.io import read
import matplotlib.pyplot as plt
import numpy as np


def bxb_angles(perovskite, bxb_scale = 1.4, supercell = (1,1,1), debug = False):
    """
    Computes the BXB angles, as well as their average. *NOTE* this utilizes the ASE cutoffs for
    what counts as a bond in order to compute these angles! This is a different approach than that used to
    consider what counts as a BX6 octahedra in perovskite.py!
    Returns:
        (average BXB angle, np array of individual BXB angles)
    """
    
    X = perovskite.X
    B = perovskite.B
    tmp_atoms = perovskite.atoms.repeat(supercell)
    ana = Analysis(tmp_atoms, cutoffs = natural_cutoffs(tmp_atoms, mult=bxb_scale))
    BXB_angs = ana.get_angles(B, X, B, unique=True)
    if debug:
        print(f"\nThere are {len(BXB_angs[0])} {B}-{X}-{B} bond angles\n")
    if len(BXB_angs[0]) > 0:
        BXB_ang_values = ana.get_values(BXB_angs)
    else:
        BXB_ang_values = None

    if perovskite.Bp is not None:
        BXBp_angs = ana.get_angles(B, X, perovskite.Bp, unique=True)
        if debug:
            print(f"\nThere are {len(BXBp_angs[0])} {B}-{X}-{perovskite.Bp} bond angles\n")
        if len(BXBp_angs[0]) > 0:
            BXBp_ang_values = ana.get_values(BXBp_angs)
        else:
            BXBp_ang_values = None
    else:
        BXBp_ang_values = None

    return(BXB_ang_values, BXBp_ang_values)

def _get_plottable_partial_rdf(perov, element_pairs, max_dist=10, npoints = None, mode = 'pyrovskite', ss_norm = False):
    """
    Computes and returns partial rdf for the element pairs in element_pairs. 
    Input:
        perov: Perovskite object
        element_pairs: List of len()=2 lists for pairs of elements to be computed. E.g. [["Pb", "I"], ["Pb", "Pb"], ["I", "I"]]
        max_dist: Max distance in angstroms for the rdf computation
        mode: 'ase', or 'pyrovksite'. ase bins the data then plots, pyrovskite extracts distances and then smooths with gaussian
                kernel. Both do the same thing, but pyrovskite method will look a bit cleaner, at the cost of being slower.
    Output:
        plottable_rdfs: Dictionary with keys of the form "PbI_x", "PbPb_x", "II_x", and "PbI_rdf", "PbPb_rdf", "II_rdf", corresponding
        to x-y axes respectively of the plottable rdf functions for the example element_pair input given above.

    Note this format is preferred for the serialization of the data in cases where many perovskite RDFs are to be combined.
    Note, if you request a max_dist incompatible with your simulation cell, a supercell will be created for you. This *may not be 
    your desired behavior* as long range periodicity may be overemphasized/artficial. Select max_dist accordingly if you want to
    avoid this.
    """

    if npoints is None:
        npoints = max_dist * 10

    # Set up supercell if needed
    cell = perov.atoms.cell
    cell_norms = [np.linalg.norm(cell[0]), np.linalg.norm(cell[1]), np.linalg.norm(cell[2])]
    supercell = [0,0,0]
    for idx, i in enumerate(cell_norms):
        scl = int(max_dist*2/i)+1
        supercell[idx] = scl
    tmp_atoms = perov.atoms.repeat(tuple(supercell))

    plottable_rdfs = {}
    if mode == 'ase':
        ana = Analysis(tmp_atoms)
        for pair in element_pairs:
            rdf_data = ana.get_rdf(max_dist, npoints, elements = pair, return_dists = True)
            plottable_rdfs["".join(pair)+"_x"] = rdf_data[0][1]
            plottable_rdfs["".join(pair)+"_rdf"] = rdf_data[0][0]

    elif mode == 'pyrovskite':

        # !!! This routine is largely taken from: https://gitlab.com/ase/ase/-/blob/master/ase/geometry/rdf.py
        # This routine does not bin the data as in the link above, but rather returns the rdf spectra for 
        # smoothing with gaussian kernel. The ase routine also appears to have a bug where it's picking up
        # certain bond distances that aren't really there. Those aren't present in this routine
        dm = tmp_atoms.get_all_distances(mic=True)
        curr_max = 0
        for pair in element_pairs:
            pair_dists = []
            i_indices = np.where(tmp_atoms.symbols == pair[0])[0]
            phi = len(i_indices)
            # is the len(i_indices part, as this changes across element_pairs.)
            for i in i_indices:
                for j in np.where(tmp_atoms.symbols == pair[1])[0]:
                    if i != j:
                        pair_dists.append(dm[i, j])

        # !!! This routine is largely taken from: https://gitlab.com/ase/ase/-/blob/master/ase/geometry/rdf.py

            this_xs, this_spectra = _gaussian_kernel_discrete_spectrum(pair_dists, smearing = .05, 
                                                                       gridpoints = max_dist * 50, v_min = 0, v_max = max_dist)
            if ss_norm:
                this_spectra *= phi

            plottable_rdfs["".join(pair)+"_x"] = this_xs
            plottable_rdfs["".join(pair)+"_rdf"] = this_spectra

            if curr_max < np.max(this_spectra):
                curr_max = np.max(this_spectra)

        for pair in element_pairs:
            plottable_rdfs["".join(pair)+"_rdf"] /= curr_max
    else:
        raise ValueError("Mode must be 'pyrovskite' or 'ase'.")
    return(plottable_rdfs)

def _gaussian_kernel_discrete_spectrum(spectrum, smearing = None, gridpoints = 200, v_min = None, v_max = None):

    """
    Smooths the input spectrum with a gaussian kernel over the supported
    region and returns the plottable information.

    Parameters
    --

    """
    spectrum = np.array(spectrum)
    if smearing is None:
        smearing = (np.max(spectrum) - np.min(spectrum))/ 50

    if v_min is None:
        v_min = np.min(spectrum) - 3*smearing
    if v_max is None:
        v_max = np.max(spectrum) + 3*smearing

    xs = np.linspace(v_min, v_max, gridpoints)
    res = np.zeros_like(xs)
    for delta in spectrum:
        res += np.exp(-(xs-delta)**2/smearing**2)

    return(xs, res)


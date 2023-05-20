from pymatgen.io.ase import AseAtomsAdaptor as ase2pmg
from pyrovskite.geometry import bxb_angles, _gaussian_kernel_discrete_spectrum, _get_plottable_partial_rdf
from pyrovskite.io import xTB_input, GPAW_input
from pymatgen.core.structure import Structure
import matplotlib, itertools, sys
import matplotlib.pyplot as plt
from ase import neighborlist
import ase, ase.io
import numpy as np

class Perovskite:
    """
    Constructor arguments:
        atoms (str or ase.atoms.Atoms): Filename or Atoms object of the 
                                        perovskite system.
        Ap (str): Name of spacer cation in 2D perovskite systems.
        A  (str): A-site cation (see list of A-site cations code will look for
                  automatically).
        B  (str): B-site cation (see list of B-site cations code will look for
                  automatically).
        Bp (str): B'-site cation in the case of Double perovskite systems.
        X  (str): X-site anion (code will expect Halogen or Oxygen).
        n  (int): Layer number in the case of 2D perovskite systems, optionally
                  use np.inf for bulk.
    """

    def __init__(self, atoms, Ap=None, A=None, B=None, Bp=None, X=None, n=None):
        if type(atoms) == str:
            self.atoms=ase.io.read(atoms)
        elif type(atoms) == ase.atoms.Atoms:
            self.atoms = atoms

        self.structure = ase2pmg.get_structure(self.atoms)
        self.Ap = Ap
        self.A  = A
        self.B  = B
        self.Bp = Bp
        self.X  = X
        self.n  = n

        # Set by self.identify_octahedra()
        self.octahedra = None
        self.bx_distances = None

        # Set by self.compute_sigma()
        self.octahedra_sigma = None
        self.sigma = None
        self.cis_angles = None
        self.trans_angles = None
        self.trans_pairs = None

        # Set by self.compute_delta()
        self.octahedra_delta = None
        self.delta = None

        # Set by self.compute_lambda()
        self.octahedra_Lambda_3 = None
        self.octahedra_Lambda_2 = None
        self.octahedra_Lambda = None
        self.Lambda_3 = None
        self.Lambda_2 = None
        self.Lambda = None

        # These are used to help determine B-, X-site ions insofar as they're not provided.
        # This can also be passed as an argument if one is using perovskites containing more 
        # exotic materials, e.g. organic X-site anions, or TM B-site cations, etc.
        self.common_B = ['Pb', 'Sn', 'Ge', 'Bi', 'In', 'Tl', 'Zn', 'Cu', 'Mn', 'Sb', 'Cd', 'Fe', 'Ag', 'Au']
        self.common_X = ['O', 'F', 'Cl', 'Br', 'I']

        unique_types = []
        for atom in self.atoms.get_chemical_symbols():
            if atom not in unique_types:
                unique_types.append(atom)

        self.atom_types = unique_types

        if self.Ap == None:
            pass
        if self.A == None:
            pass
        if self.B == None:
            B_candidates = list(set(self.common_B) & set(self.atom_types))
            print(f"\nNo B-cation set, candidate B-cations determined from structure:\n{B_candidates}\n")

            if len(B_candidates) == 1:
                self.B = B_candidates[0]
            elif len(B_candidates) == 2:
                print(f"Two B-site candidates found, possible double perovskite system\nSetting B={B_candidates[0]} and Bp={B_candidates[1]}")
                self.B = B_candidates[0]
                self.Bp = B_candidates[1]
            else: # Figure out how to deal with weird structures later.
                ...
                sys.exit(1)

        if self.X == None:
            # The X anions are more problematic because they can pop up in spacing cations.
            # Particularly O, F.
            self.common_X = ['O', 'F', 'Cl', 'Br', 'I']
            X_candidates = list(set(self.common_X) & set(self.atom_types))
            print(f"\nNo X-anion set, candidate X-anions determined from structure:\n{X_candidates}\n")
            if len(X_candidates) == 1:
                self.X = X_candidates[0]

            # Here we look at fragments, and try to find those most commonly bound with a B-site cation
            # The logic here is this: 
            # For spacers containing, e.g. F, O, their nearest neighbors should be C, N, O, H.
            # For F, O that are true X-site cations, their nearest neighbors should be B-site cations.
            else:
                X_candidates = list(set(self.common_X) & set(self.atom_types))
                org_candidates = ['C', 'H', 'N', 'O']

                # Convert Atoms(ASE) -> Structure(pymatgen)
                # Distance matrix with nonzero diagonal.
                dist_mat_nz_diag = self.structure.distance_matrix.copy()
                np.fill_diagonal(dist_mat_nz_diag, 100) # Allows min(row_of_dist_matrix) to be used.

                # Setup the counts of nearest neighbors for all the X candidates.
                X_nn = {X:{'inorg':0, 'org':0} for X in X_candidates}

                # Update nearest neighbor type (organic vs inorganic) for the X candidate atoms.
                for idx, site in enumerate(self.structure.sites):
                    if str(site.specie) in X_candidates:
                        nn_dist = min(dist_mat_nz_diag[idx])
                        nn_idx  = np.where(dist_mat_nz_diag[idx] == nn_dist)
                        nn_spec = str(self.structure.sites[nn_idx[0][0]].specie)
                        if nn_spec in org_candidates:
                            X_nn[str(site.specie)]['org'] +=1
                        elif nn_spec in self.common_B or nn_spec in X_candidates:
                            X_nn[str(site.specie)]['inorg'] +=1


                X_max = None
                inorg_max = 0
                # Use local chemical environment -> X_max is suggested anion.
                for key in X_nn.keys():
                    if X_nn[key]['inorg'] > inorg_max:
                        X_max = key
                        inorg_max = X_nn[key]['inorg']
                




                # If local chem env. determines X-site anion.
                # If not, we need a slightly more robut (but more time)
                # consuming approach to determining the X-anion
                if X_max is not None:
                    print(f"{X_max} detected as the X-anion. [From local env.]")
                    self.X = X_max
                else:
                    # Approach 2:
                    # For each X_c in X_candidate:
                    #   Loop through all X_c sites.
                    #   Find closest B/Bp site.
                    # Average these by X_candidate
                    # Select candidate with lowest avg.
                    # nearest B-/Bp-cation neighbor distance.
                    X_nb = {X:[] for X in X_candidates}
                    for idx, site in enumerate(self.structure.sites):
                        if str(site.specie) in X_candidates:
                            curr_min = 500
                            # Pull closest B/Bp-cation for each instance of X_candidate.
                            for idx2, site2 in enumerate(self.structure.sites):
                                if str(site2.specie) == self.B or str(site2.specie) == self.Bp:
                                    curr_distance = dist_mat_nz_diag[idx, idx2]
                                    if curr_distance < curr_min:
                                        curr_min = curr_distance
                            X_nb[str(site.specie)].append(curr_min)

                    X_id = list(X_nb.keys())
                    X_av_list = [sum(nb_list)/len(nb_list) for nb_list in X_nb.values()]

                    tmp_idx = X_av_list.index(min(X_av_list))
                    self.X = X_id[tmp_idx]
                    print("\nSimple local env method failed.")
                    print(f"{self.X} determined as X-anion from more robust method. Check this is right.\n")

                # STOICHIOMETRY NEEDS MORE WORK... 2D Perovskites have X_3n+1, bulk have X_3n, so we
                # Need to be able to differentiate between the two.. 
                # # Use stoichiometry -> stoich_X is suggested anion
                # B_ct = sum([0 if x.symbol != self.B and x.symbol != self.Bp else 1 for x in self.atoms])
                # X_cts = np.zeros(len(X_candidates))
                # for Xc_idx, Xc in enumerate(X_candidates):
                #     X_cts[Xc_idx] = sum([0 if x.symbol != Xc else 1 for x in self.atoms])
                # stoich_X = X_candidates[np.argmin(X_cts % 3)]
                # print(f"{stoich_X} detected as the X-anion. [From stoichiometry, local env. failed]")
                    

    def identify_octahedra(self, return_distances = False):
        """
        Args:
            return_distance (bool): Whether or not to return the B-X distance of the sites.
        Returns:
            octahedra (list): List of lists of the format:
                 [[ Bx,  By,  Bz],
                  [X1x, X2x, X1z],
                  [X2x, X2y, X2z],
                  [X3x, X3y, X3z],
                  [X4x, X4x, X4z],
                  [X5x, X5y, X5z],
                  [X6x, X6y, X6z]]
               where the X1 is closests to B and X6 is farthest.
               E.g. to access the 5th X-atom's z-position from the 2nd octahedra, you'd do:
               octahedra[1][4][2]

            distances (list): List of the format:
                [[B-X1 dist, B-X2 dist, ... , B-X6 dist], ... for each octahedra]
        """

        B = self.B
        Bp = self.Bp
        X = self.X

        octahedra = []
        octa_indices = []
        distances = []

        # Loop through all sites
        # This whole codeblock may be removed soon.
        pmg_site_list = []
        for i in range(len(self.atoms)):
            pmg_site_list.append([self.structure[i].x, self.structure[i].y, self.structure[i].z])
        pmg_site_list = np.array(pmg_site_list)

        for idx, site in enumerate(self.structure.sites):

            # Find the B-cations as the anchoring point for octahedra.
            if str(site.specie) == B or str(site.specie) == Bp:
                octahedra.append([site.coords]) # Append coordinates of B-cation
                octa_indices.append([idx])      # Append index of B-cation
                oct_ct = len(octahedra) 
                neighbors = self.structure.get_neighbors(site, 10.0)
                X_neighbors = []
                for neighbor in neighbors:
                    # Neighbor can be accessed as 4-tuple with attributes:
                    # (site, distance, index, image)
                    if str(neighbor.specie) == X:
                        X_neighbors.append([neighbor[0], neighbor[1], neighbor[3]])

                # Sorts and then slices, such that only 6 closest X-atoms are present.
                sorted_X_neighbors = sorted(X_neighbors, key=lambda x: x[1])[:6]
                distances.append([s[1] for s in sorted_X_neighbors])
                octahedra[oct_ct-1].extend([xi[0].coords for xi in sorted_X_neighbors])

                # May be deprecated
#                print(pmg_site_list)
#                #print([np.argmin()self.structure.sites.index([xsite[0].x, xsite[0].y, xsite[0].z]) for xsite in sorted_X_neighbors])
#                print([np.argmin(np.linalg.norm(pmg_site_list - np.array([xsite[0].x, xsite[0].y, xsite[0].z]), axis = 1))  \
#                        for xsite in sorted_X_neighbors])
#                print(np.linalg.norm(pmg_site_list - np.array([sorted_X_neighbors[0][0].x, sorted_X_neighbors[0][0].y, sorted_X_neighbors[0][0].z]), axis = 1))
        self.octahedra = np.array(octahedra)
        self.bx_distances = np.array(distances)

        if return_distances:
            return(np.array(octahedra), np.array(distances))
        else:
            return(np.array(octahedra))

    def plot_vertices(self, show = True, save = False, filename = None):
        if save and filename is None:
            filename = f"{self.B}_{self.X}_octahedra.png"

        if self.octahedra is None:
            self.identify_octahedra()
        if self.trans_pairs is None:
            self.compute_sigma()

        max_bx = np.max(self.bx_distances) ** 2 * 3.14 * 50
        fig = plt.figure(12)
        ax = fig.add_subplot(projection='3d')

        c_dicts = []
        for idx, octa in enumerate(self.octahedra):
            c_dicts.append({})
            colors = ['red', 'green', 'blue']
            for pair in range(3):
                c_dicts[idx][self.trans_pairs[idx,pair,0]] = colors[pair]
                c_dicts[idx][self.trans_pairs[idx,pair,1]] = colors[pair]

            ax.scatter(octa[0,0], octa[0,1], octa[0,2], marker = '8', s=max_bx*.35)
            for x_idx, pt in enumerate(octa[1:]):
                this_color = c_dicts[idx].get(x_idx, 'black')
                ax.plot([octa[0,0], pt[0]], [octa[0,1], pt[1]], [octa[0,2], pt[2]], color=this_color)

        xs = self.octahedra[:,1:,0]
        ys = self.octahedra[:,1:,1]
        zs = self.octahedra[:,1:,2]
        ax.scatter(xs, ys, zs, marker='o', s=max_bx*.15)
        plt.xlim()

        if save:
            plt.savefig(filename, dpi=500)
        if show:
            plt.show()


    def _get_angle(self, a, b):
        """ Simple trig to get X-B-X angle"""
        return np.arccos((np.dot(a,b)/(np.linalg.norm(a)*np.linalg.norm(b)))) * 180/np.pi

    def compute_delta(self, return_type = "delta"):
        """
        Computes: 1/6 Sigma_{BX} | ((d_BX - d_m)/d_m)^2 |
        Returns:
            return_type = "delta" -> Average delta for all octahedra
            return_type = "octahedra_delta" -> Delta for each octahedra
            return_type = "both" -> (delta, octahedra_delta)
            Also sets both self.delta and self.octahedra_delta
        """
        if self.octahedra is None:
            self.identify_octahedra()

        octahedra_delta = np.zeros(self.bx_distances.shape[0])
        for idx, octa_distances in enumerate(self.bx_distances):
            this_octa_avg = np.mean(octa_distances)
            octahedra_delta[idx] = np.mean(((octa_distances - this_octa_avg)/this_octa_avg)**2) 

        self.delta = np.mean(octahedra_delta)
        self.octahedra_delta = octahedra_delta

        if return_type == "delta":
            return(self.delta)
        elif return_type == "octahedra_delta":
            return(self.octahedra_delta)
        elif return_type == "both":
            return(self.delta, self.octahedra_delta)
        else:
            print("Invalid return_type, enter delta, octahedra_delta, or both. See help(Perovskite.compute_delta()) for info")

    def compute_sigma(self, return_type = "sigma"):
        """
        Computes: 1/12 * Sigma_{cis} | phi_{cis} - 90^o |
        Sigma is an octahedra specific parameter, so much of this function is hardcoded, and it
        is not applicable to other coord complex polyhedra

        Returns:
            return_type = "sigma" -> Average sigma for all octahedra
            return_type = "octahedra_sigma" -> Sigma for each octahedra
            return_type = "both" -> (sigma, octahedra_sigma)
            Also sets both self.sigma and self.octahedra_sigma
        """
        # sigma is specific to octahedra, so a lot of things are hardcoded.
        if self.octahedra is None:
            self.identify_octahedra()

        # Another comment on pair_list, there are 12-cis angles, and 3-trans angles.
        # Instead of having a complex geometric evaluation of which the cis-trans angles are,
        # we simply remove the 3 angles closest to 180 degrees.
        pair_list = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (1, 3), (1, 4), 
                     (1, 5), (2, 3), (2, 4), (2, 5), (3, 4), (3, 5), (4, 5)]

        cis_angle_list = np.zeros((self.octahedra.shape[0], 12))
        trans_angle_list = np.zeros((self.octahedra.shape[0], 3))
        all_angles = np.zeros((self.octahedra.shape[0], 15))
        all_angles_pairs = np.zeros((self.octahedra.shape[0], 15, 2), dtype='int')
        #                               Otctahedra, trans_idx, X-idxs
        all_trans_pairs = np.zeros((self.octahedra.shape[0], 3, 2), dtype='int')
        octahedra_sigma = np.zeros(self.octahedra.shape[0])
        for idx, octa in enumerate(self.octahedra):
            # Get relative positions with B-cation at the origin.
            B_pos = octa[0] - octa[0]
            X_pos = octa[1:] - octa[0]
            # Assert that bx bond distances didn't change upon shifting to CoM
            assert np.allclose(np.linalg.norm(X_pos, axis = 1),self.bx_distances[idx]) , "BX Bond distances before and after shift to relative coordinates are not the same!"
            
            for idx2, pair in enumerate(pair_list):
                this_angle = self._get_angle(X_pos[pair[0]], X_pos[pair[1]])
                all_angles[idx, idx2] = this_angle
                all_angles_pairs[idx,idx2,0] = pair[0]
                all_angles_pairs[idx,idx2,1] = pair[1]

            # This gets us the 3 angles closest to 180 (our trans-angles)
            closest_to_180 = np.sort(np.abs(all_angles[idx]-180.0) )[:3]
            close_plus  = 180.0 + closest_to_180
            close_minus = 180.0 - closest_to_180

            cis_idx = 0
            trans_idx = 0
            # Loop through all_angles, and remove the three closest to 180.
            # This fills out cis_angle_list with our cis angles.
            for angles_idx, x in enumerate(all_angles[idx]):
                if True in np.isclose(close_plus, x) or True in np.isclose(close_minus, x):
                    trans_angle_list[idx, trans_idx] = x
                    all_trans_pairs[idx, trans_idx, 0] = all_angles_pairs[idx, angles_idx,0]
                    all_trans_pairs[idx, trans_idx, 1] = all_angles_pairs[idx, angles_idx,1]
                    trans_idx += 1
                else:
                    cis_angle_list[idx, cis_idx] = x
                    cis_idx += 1
                assert cis_idx <= 12 , "too many elements being added to cis_angle_list"

            octahedra_sigma[idx] = np.mean(np.abs(cis_angle_list[idx] - 90.0))

        self.octahedra_sigma = octahedra_sigma
        self.sigma = np.mean(octahedra_sigma)
        self.cis_angles = cis_angle_list
        self.trans_angles = trans_angle_list
        self.trans_pairs = all_trans_pairs

        if return_type == "sigma":
            return(self.sigma)
        elif return_type == "octahedra_sigma":
            return(self.octahedra_sigma)
        elif return_type == "both":
            return(self.sigma, self.octahedra_sigma)
        else:
            print("Invalid return_type, enter sigma, octahedra_sigma, or both. See help(Perovskite.compute_sigma()) for info")

    def compute_lambda(self, visualize=False, return_type = "lambda", scaled = True, ratio = False):

        # Need the trans atom indices from compute_sigma.
        if self.trans_pairs is None:
            self.compute_sigma()

        if visualize:
            fig = plt.figure(12)
            ax = fig.add_subplot(projection='3d')

        # Do lambda computation for all octahedra.
        Lambda_3_list = np.zeros(self.octahedra.shape[0])
        Lambda_2_list = np.zeros(self.octahedra.shape[0])
        Lambda_list = np.zeros(self.octahedra.shape[0])
        for idx, octa in enumerate(self.octahedra):
            B = octa[0]         # B-site cation
            Xs = octa[1:]       # X-site anions
            tilde_P_list = []   # Holder for midpoints of the X-X connections.

            res = np.zeros_like(Xs[0])  # For computation of tilde{P} from tilde_P_list
            e_basis = np.zeros((3,3))   # For the 'basis' of the octahedra
            print(e_basis)
            for e_idx, pair in enumerate(self.trans_pairs[idx]):
                
                # Here we store the basis vectors for the octahedra.
                e_basis[e_idx] = ((Xs[pair[0]] - Xs[pair[1]]) / 
                np.linalg.norm(Xs[pair[0]] - Xs[pair[1]]))
                
                # Compute bisection of the X-X connection
                mid_point = 0.5*(Xs[pair[0]]+Xs[pair[1]])
                tilde_P_list.append(mid_point)
                res += mid_point

            # Compute P_tilde, the 'natural center' of a distorted octahedra.
            # For a perfect octahedra, this directly coincides with the B-site cation position.
            P_tilde = res / 3

            # Displacement vector of the B-cation from its 'natural center'
            D = B-P_tilde
            D_i = np.array([abs(np.dot(D,e_basis[0])), abs(np.dot(D, e_basis[1])),
                            abs(np.dot(D, e_basis[2]))])

            s=f"""e_basis not unit length\n|D|: {np.linalg.norm(D)}\n|D_i|: {np.linalg.norm(D_i)}
            \n|e1| = {np.linalg.norm(e_basis[0])}\n|e2| = {np.linalg.norm(e_basis[1])}\n|e3| = {np.linalg.norm(e_basis[2])}"""
            #assert abs(np.linalg.norm(D) -  np.linalg.norm(D_i)) < 1e-5 ,  s

            lambda_ij = np.zeros(3)

            if scaled:
                lambda_ij[0] = min(D_i[0]/D_i[1], D_i[1]/D_i[0]) * (D_i[0]+D_i[1])
                lambda_ij[1] = min(D_i[0]/D_i[2], D_i[2]/D_i[0]) * (D_i[0]+D_i[2])
                lambda_ij[2] = min(D_i[1]/D_i[2], D_i[2]/D_i[1]) * (D_i[2]+D_i[1])
            else:
                lambda_ij[0] = min(D_i[0]/D_i[1], D_i[1]/D_i[0]) 
                lambda_ij[1] = min(D_i[0]/D_i[2], D_i[2]/D_i[0])
                lambda_ij[2] = min(D_i[1]/D_i[2], D_i[2]/D_i[1])

            # Key parameters.
            Lambda_3_list[idx] = lambda_ij[0]*lambda_ij[1]*lambda_ij[2]
            Lambda_2_list[idx] = np.max(lambda_ij)
            Lambda_list[idx] = Lambda_3_list[idx]/Lambda_2_list[idx]

            # Visualization, might remove later.
            if visualize:
                # Ensures reasonably sized markers.
                max_bx = np.max(self.bx_distances) ** 2 * 3.14 * 50
                # Plots the B-cation site.
                ax.scatter(B[0], B[1], B[2], marker = '8', s=max_bx*.05, label='B')
                # Plots the 'natural center' P_tilde
                ax.scatter(P_tilde[0], P_tilde[1], P_tilde[2], marker = '^',
                           s=max_bx*.05, label='$\\tilde{P}$')
                # Plots the bisection points of the three X-X connections.
                for P in tilde_P_list:
                    ax.scatter(P[0], P[1], P[2], marker = 'x', s=max_bx*.05)

                # Plots the 'basis vectors' of the octahedra.
                ax.plot([P_tilde[0], (P_tilde+e_basis[0])[0]],
                        [P_tilde[1], (P_tilde+e_basis[0])[1]],
                        [P_tilde[2], (P_tilde+e_basis[0])[2]], color='black')
                ax.plot([P_tilde[0], (P_tilde+e_basis[1])[0]],
                        [P_tilde[1], (P_tilde+e_basis[1])[1]],
                        [P_tilde[2], (P_tilde+e_basis[1])[2]], color='black')
                ax.plot([P_tilde[0], (P_tilde+e_basis[2])[0]],
                        [P_tilde[1], (P_tilde+e_basis[2])[1]],
                        [P_tilde[2], (P_tilde+e_basis[2])[2]], color='black')

        self.octahedra_Lambda_3 = Lambda_3_list
        self.octahedra_Lambda_2 = Lambda_2_list
        self.octahedra_Lambda = Lambda_list
        self.Lambda_3 = np.average(Lambda_3_list)
        self.Lambda_2 = np.average(Lambda_2_list)
        self.Lambda = np.average(Lambda_list)

        if visualize:
            plt.title("x = Midpoint of trans X-X bond\no = B-cation\n^ = $\\tilde{P}$")
            plt.show()

        if ratio:
            if return_type == "lambda":
                return(self.Lambda_3, self.Lambda_2, self.Lambda)
            elif return_type == "octahedra_lambda":
                return(self.octahedra_Lambda_3, self.octahedra_Lambda_2, self.octahedra_Lambda)
            elif return_type == "both":
                print("compute_lambda return both signature is:\n octa_L3, octa_L2, octa_L, L3, L2, L")
                return(self.octahedra_Lambda_3, self.octahedra_Lambda_2, self.octahedra_Lambda,
                       self.Lambda_3, self.Lambda_2, self.Lambda)
            else:
                print("Invalid return type, select lambda, octahedra_lambda, or both.")
        else:
            if return_type == "lambda":
                return(self.Lambda_3, self.Lambda_2)
            elif return_type == "octahedra_lambda":
                return(self.octahedra_Lambda_3, self.octahedra_Lambda_2)
            elif return_type == "both":
                print("compute_lambda return both signature is:\n octa_L3, octa_L2, octa_L, L3, L2, L")
                return(self.octahedra_Lambda_3, self.octahedra_Lambda_2,
                       self.Lambda_3, self.Lambda_2)
            else:
                print("Invalid return type, select lambda, octahedra_lambda, or both.")

    def plot_angles(self, smearing = 1, gridpoints = 300, fignum = 123,
                    show = True, bxb_scale = 1.5, supercell = (1,1,1), debug = False):

        if self.octahedra is None:
            self.identify_octahedra()
        if self.cis_angles is None or self.trans_angles is None:
            self.compute_sigma()

        # Compute the BXB angles (Not needed for any octahedra distortion parameters)
        BXB_angs, BXBp_angs = bxb_angles(self, bxb_scale = bxb_scale,
                                         supercell = supercell, debug = debug)

        # Plotting
        plt.figure(fignum, figsize = (12,6))
        font = {'size'   : 18}
        matplotlib.rc('font', **font)


        plt.title("Perovskite Bond Angle Distributions", font = font)
        if BXB_angs is not None:
            BXB_angs = BXB_angs[0]
            # Signature:
            # _gaussian_kernel_discrete_spectrum(spectrum, smearing = None, gridpoints = 200):
            xs, smeared_BXB_angs = _gaussian_kernel_discrete_spectrum(BXB_angs, smearing = smearing, gridpoints = gridpoints)
            smeared_BXB_angs *= 1/(supercell[0]*supercell[1]*supercell[2])
            plt.plot(xs, smeared_BXB_angs, label = "B-X-B", color = 'purple', linewidth = 2)

        if BXBp_angs is not None:
            BXBp_angs = BXBp_angs[0]
            xps, smeared_BXBp_angs = _gaussian_kernel_discrete_spectrum(BXBp_angs, smearing = smearing, gridpoints = gridpoints)
            smeared_BXBp_angs *= 1/(supercell[0]*supercell[1]*supercell[2])
            plt.plot(xps, smeared_BXBp_angs, label = "B-X-Bp",
                     color = 'black', linewidth = 2)

        xts, smeared_trans_angs =  _gaussian_kernel_discrete_spectrum(np.ravel(self.trans_angles), 
                                                                      smearing = smearing, gridpoints = gridpoints)
        xcs, smeared_cis_angs =  _gaussian_kernel_discrete_spectrum(np.ravel(self.cis_angles),
                                                                    smearing = smearing, gridpoints = gridpoints)
        if self.Bp is not None:
            plt.plot(xts, smeared_trans_angs, label = "Trans X-B/Bp-X",
                     color = 'red', linewidth = 2)
            plt.plot(xcs, smeared_cis_angs, label = "Cis X-B/Bp-X", color = 'blue',
                     linewidth = 2)
        else:
            plt.plot(xts, smeared_trans_angs, label = "Trans X-B-X", color = 'red',
                     linewidth = 2)
            plt.plot(xcs, smeared_cis_angs, label = "Cis X-B-X", color = 'blue',
                     linewidth = 2)

        plt.ylim(ymin=.01)
        plt.xlabel("Bond Angle [°]")
        plt.ylabel("Intensity [arb. u.]")
        plt.legend()
        if show:
            plt.show()


    def plot_distances(self, smearing = .02, gridpoints = 300, fignum = 12,
                       show = True):

        # Ensuring data is here for plotting
        if self.octahedra is None:
            self.identify_octahedra()
        if self.bx_distances is None:
            self.compute_delta()

        # Plotting
        plt.figure(fignum, figsize = (12,6))
        font = {'size'   : 18}
        matplotlib.rc('font', **font)

        xs, smeared_BX_dists = _gaussian_kernel_discrete_spectrum(np.ravel(self.bx_distances),
                                                                  smearing = smearing,
                                                                  gridpoints = gridpoints)

        if self.Bp is None:
            plt.title(f"{self.B}-{self.X} Bond Distances")
        else:
            plt.title(f"{self.B}-{self.X}/{self.Bp}-{self.X} Bond Distances")

        plt.plot(xs, smeared_BX_dists, color = 'black', linewidth = 2)
        plt.ylim(ymin=.01)
        plt.xlabel("Distance [Å]")
        plt.ylabel("Intensity [arb. u.]")
        if show:
            plt.show()

    def plot_rdf(self, max_dist = 10, ss_norm = False, mode = 'pyrovskite',
                 fignum = 22, show = True):
        """
        Input:
            max_dist (float/int): Max distance in angstrom for the RDF 
                                  computation
            ss_norm (bool): Whether or not to use normalization which is propto
                            the # of atoms in the computation, either way max
                            peak is scaled to 1.
        Output:
            Nothing, only plot.

        Note, if you request a max_dist incompatible with your simulation cell,
        a supercell will be created for you. This *may not be your desired 
        behavior* as long range periodicity may be overemphasized/artficial.
        Select max_dist accordingly if you want toavoid this.
        """
        if self.octahedra is None:
            self.identify_octahedra()

        if self.Bp is None:
            element_pairs = [[self.B, self.X],
                             [self.B, self.B],
                             [self.X, self.X]]
        else:
            element_pairs = [[self.B, self.X],
                             [self.B, self.B],
                             [self.X, self.X],
                             [self.Bp, self.X],
                             [self.Bp, self.Bp],
                             [self,Bp, self.B]]

        plottable_rdfs = _get_plottable_partial_rdf(self, element_pairs,
                                                    max_dist = max_dist,
                                                    mode = mode,
                                                    ss_norm=ss_norm)

        plt.figure(fignum, figsize = (12,6))
        font = {'size'   : 18}
        matplotlib.rc('font', **font)

        for key in element_pairs:
            plt.plot(plottable_rdfs["".join(key)+"_x"],
                     plottable_rdfs["".join(key)+"_rdf"],
                     label = key[0]+"-"+key[1])
        plt.title("Partial Radial Distribution Function")
        plt.xlim(0, max_dist)
        plt.ylim(ymin = 0.01)
        plt.legend()
        if show:
            plt.show()

    def write_xTB(self, calc_type = "vcopt", prefix = "xtb", cif_name = None,
                  atoms = None, xtb_path = None, d3_path = None):
        """ 

        This function is a thin wrapper around pyrovskite.input_generator's
        xTB_input function. See help(pyrovskite.input_generator.xTB_input)
        for the full documentation.

        """

        if atoms is None:
            atoms = self.atoms

        xTB_input(calc_type = calc_type, prefix = prefix, cif_name = cif_name,
                  atoms = atoms, xtb_path = xtb_path, d3_path = d3_path)

    def write_GPAW(self, calc_type = "fcopt", prefix = "gpaw", cif_name = None,
                   atoms = None, opt_algo = "FIRE", fmax = .02):
        """ 

        This function is a thin wrapper around pyrovskite.input_generator's
        GPAW_input function. See help(pyrovskite.input_generator.GPAW_input)
        for the full documentation.

        """

        if atoms is None:
            atoms = self.atoms
        
        GPAW_input(calc_type = "fcopt", prefix = "gpaw", cif_name = None,
              atoms = atoms, opt_algo = "FIRE", fmax = .02)


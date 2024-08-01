
def get_ionic_radius(site, name):
    site_dict = ionic_radii.get(site)
    val = site_dict.get(name, None)

    if val is not None:
        return val
    else:
        print_available_ionic_radii()
        print("\nIon not available, available ions printed above")
        raise ValueError("Ion requested does not have info in the code. See the above message")


def print_available_ionic_radii(mode = "all"):
    """
    Prints all available ionic radii.
    """

    print("===================================NOTE==========================================")
    print("|  All B-site cations are assumed to be divalent unless otherwise specified.    |")
    print("|                        All X-site halogens are -1.                            |")
    print("|                        All A-site cations are +1.                             |")
    print("| References for all listed radii can be found in pyrovskite/pyrovskite/data.py |")
    print("===================================NOTE==========================================\n")


    if mode == "all":
        print("============")
        print("|A cations.|")
        print("============")
        for a in ionic_radii["A"].keys():
            print(f"{a}:\tAlternative names: {A_cation_aliases[a]}")
        print("\n============")
        print("|B cations.|")
        print("============")
        print(list(ionic_radii["B"].keys()))
        print("\n==========")
        print("|X anions.|")
        print("===========")
        print(list(ionic_radii["X"].keys()))


# References for the listed data:
# [1] https://doi.org/10.1021/acs.chemrev.8b00539
# [2] https://doi.org/10.1021/acs.chemmater.9b05273
# [3] https://doi.org/10.1107/S0567739476001551

ionic_radii = {
    "A":{

        # This is completely arbitrary choice of Cs for the perovskites
        # that don't contain an A-site cation (e.g. n = 1 2DPKs).
        "NA"   : 1.88, # See Note above

        "NH4"  : 1.46, # [1]
        "MA"   : 2.17, # [1]
        "FA"   : 2.53, # [1]
        "HZA"  : 2.17, # [1]
        "AZ"   : 2.50, # [1]
        "HXA"  : 2.16, # [1]
        "IMA"  : 2.58, # [1]
        "EA"   : 2.74, # [1]
        "DMA"  : 2.72, # [1]
        "GA"   : 2.78, # [1]
        "TMA"  : 2.92, # [1]
        "TA"   : 3.20, # [1]
        "3-PYR": 2.72, # [1]
        "TPY"  : 3.33, # [1]
        "K"    : 1.64, # [1]
        "Rb"   : 1.72, # [1]
        "Cs"   : 1.88, # [1]
        "MHy"  : 2.64, # [2]
    },
    "B":{
        "Pb"   : 1.19, # [1]
        "Sn"   : 1.10, # [1]
        "Ge"   : 0.73, # [1]
        "Mg"   : 0.72, # [1]
        "Ca"   : 1.00, # [1]
        "Sr"   : 1.18, # [1]
        "Ba"   : 1.35, # [1]
        "Cu"   : 0.73, # [1]
        "Fe"   : 0.78, # [1]
        "Pd"   : 0.86, # [1]
        "Eu"   : 1.17, # [1]
        "Bi"   : 1.03, # [1]
        "Sb3+" : 0.76, # [1]
        "Co"   : 0.79, # [3] 
        "Hg"   : 1.16, # [3]
        "Zn"   : 0.88, # [3]
        "Cd"   : 1.09, # [3]
    },
    "X":{
        "F"    : 1.29, # [1]
        "Cl"   : 1.81, # [1]
        "Br"   : 1.96, # [1]
        "I"    : 2.20, # [1]
    },
}

A_cation_aliases = {
    "NH4"  : ["ammonium", "ammonium cation"],
    "MA"   : ["methylammonium", "[CH3NH3]+", "CH3NH3"],
    "FA"   : ["formamidinium", "[CH(NH2)2]+", "CH(NH2)2"],
    "HZA"  : ["hydrazinium", "[NH3NH2]+", "NH3NH2"],
    "AZ"   : ["azetidinium", "[(CH2)3NH2]+", "(CH2)3NH2"],
    "HXA"  : ["hydroxylammonium", "[NH3OH]+", "NH3OH"],
    "IMA"  : ["imidazolium", "[C3N2H5]+", "[C3N2H5]"],
    "EA"   : ["ethylammonium", "[(CH3CH2)NH3]+", "(CH3CH2)NH3"],
    "DMA"  : ["dimethylammonium", "[(CH3)2NH2]+", "(CH3)2NH2"],
    "GA"   : ["guanidinium", "[(NH2)3C]+", "(NH2)3C"],
    "TMA"  : ["tetramethylammonium", "[(CH3)4N]+", "(CH3)4N"],
    "TA"   : ["thiazolium", "[C3H4NS]+", "C3H4NS"],
    "3-PYR": ["3-pyrrolinium", "[NC4H8]+", "NC4H8"],
    "TPY"  : ["tropylium", "[C7H7]+", "C7H7"],
    "MHy"  : ["methylhydrazinium", "[CH7N2]+", "CH7N2"],
    "K"    : ["K+", "potassium"],
    "Cs"   : ["Cs+", "cesium"],
    "Rb"   : ["Rb+", "rubidium"],
}

import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.analysis import distances
import sys
import pandas as pd
from matplotlib import pyplot as plt

START_TIME = 0
STOP_TIME = 50

# the protein key is used to select residues in all monomers
# the protein is parsed into different monomers automatically
PROTEIN_KEY = "resid 1-50"
PROTEIN_FITTING_CRITERIA = "name CA"

# if more ligands are needed, add them using the reffered notation
# specific atom IDs can selected for the ligands 
LIGANDS_KEYS = ["resname SAM", "resname NAS"]
LIGAND_FITTING_CRITERIA = "not name H*"

#atom pair indices
ATOM_PAIRS = [(49, 50), (99, 100), (101, 103)]

FUNCTIONS = ["Distance between atoms"]
#FUNCTIONS = ["RMSD fitted ligands", "RMSD fitted protein"]

"""
Choose functions from function selection pool. Functions are applied from left to right
function options:

RMSD fitted ligands - calculates the RMSD over all frames for each inputed ligand,
fitted for each ligand

RMSD fitted protein - calculates the RMSD over all frames for the protein and 
the ligands, fitted for the inputted protein backbone

Distance between atoms - calculates the distance between 2 atoms over all frames
-> for more that an atom pair, input as a list of tuples: [(A1, A2), (A3, A4)]
"""

def analyse_rgyr(u, atom_list, max):
    """
    Computes the gyration radius for each corresponding frame and simulation time. The outputted values are dependent on the inputed step
    """
    u.trajectory[0].frame
    reduced_step = range(0,len(u.trajectory)+1)
    rgyr_data = []
    for ts in reduced_step:
        frame = ts.frame
        rgyr_data.append((u.trajectory.time, atom_list.radius_of_gyration()))
        if frame == max:
            break
    u.trajectory[0]

    export_name = "rgyr_analysis"
    rgyt_df = pd.DataFrame(rgyr_data, collumns = ["Frame", "Time", "Rgyr"])
    export_csv_plot(export_name, rgyr_data)
    return "Radius of gyration analysis sucessfully generated. \n"


def parse_protein(atom_list):
    """
    Creates a dict divided by complexes and adds a different monomer to each. For instance, a dimer dict \
    will have 2 complexes and a tetramer dict will have 4 complexes.
    
    :param atom_list: An AtomGroup object containing atoms of a protein.
    """
    complexes = {}
    current_complex = 1
    current_resid_number = -1
    current_index = -1
    for i in range(len(atom_list)):
        atom_info = atom_list[i]
        current_index = int(atom_info.index)
        resid_number = atom_info.resid
        #update the chain number if a new chain has been found
        if resid_number < current_resid_number:
            current_complex += 1
        current_resid_number = resid_number
        #add a chain to the dict if it does not exist
        if f"{current_complex}" not in complexes:
            complexes[f"{current_complex}"] = {}
        if "protein_chain" not in complexes[f"{current_complex}"]:
            complexes[f"{current_complex}"] ["protein_chain"] = [current_index]
        else:
            complexes[f"{current_complex}"] ["protein_chain"].append(current_index)

    return complexes

def add_ligands_to_complexes(u, complexes_dict, ligand_key_list):
    """
    Adds the ligands of each chain complex in the system to a pre-existing complexes-dict
    
    :param u: universe
    :param complexes_dict: a complexes dictionaty that contains protein chains for each chain complex
    :param first_lai: the first atom id of the first appearing ligand in the topology file
    """
    final_dict = complexes_dict.copy()

    for ligand_key in ligand_key_list:
        ligand_atom_group = u.select_atoms(ligand_key)
        current_resid = ligand_atom_group[0].resid
        ligand_name = ligand_atom_group[0].resname
        current_complex = 1

        for atom_info in ligand_atom_group:
            current_index = int(atom_info.index)

            if atom_info.resid > current_resid:
                current_complex += 1
                current_resid = atom_info.resid

            if "ligands" not in final_dict[f"{current_complex}"]:
                final_dict[f"{current_complex}"]["ligands"] = {}
                final_dict[f"{current_complex}"]["ligands"][ligand_name] = [current_index]

            elif ligand_name not in final_dict[f"{current_complex}"]["ligands"]:
                final_dict[f"{current_complex}"]["ligands"][ligand_name] = [current_index]

            else:
                final_dict[f"{current_complex}"]["ligands"][ligand_name].append(current_index)

        current_complex = 1
        
    return final_dict

def compute_RMSD_fit_protein(u, full_dict, complex_id, ligand_names):
    """
    Computes the RMSD of a inputed complex with fitting on the complex protein_chain, from a universe and a dictionary of all atoms in the system
    from MDAnalysis.analysis import rms
    :param u: Universe
    :param full_dict: disctionary containing atom ids, separated into complexes. Each complex contains \
    the specific monomer and ligands atoms.
    :param complex_n: the identifier of a complex
    """
    export_name = f"RMSD_Protein_fitted_plus_{'_'.join(ligand_names)}"
    complex_atom_group, complex_keys_list = get_complex_data(u, full_dict, complex_id)

    r = rms.RMSD(complex_atom_group, complex_atom_group, select = PROTEIN_FITTING_CRITERIA, groupselections = complex_keys_list)
    r.run()

    run_df = pd.DataFrame(r.results.rmsd, columns = ["Frame", "Time(ns)", "Backbone", "Alpha Carbons"] +  ligand_names)
    
    export_csv_plot(export_name, run_df)

    return "Protein fitted RMSD successfuly generated. \n"


def compute_RMSD_fit_ligand(u, full_dict, complex_id, ligand_names ):
    """
    Computes the RMSD of a inputed complex with fitting on the ligand, from a universe and a dictionary of all atoms in the system
    
    :param u: Universe
    :param full_dict: disctionary containing atom ids, separated into complexes. Each complex contains \
    the specific monomer and ligands atoms.
    :param complex_n: the identifier of a complex
    """

    for ligand in ligand_names:
        export_name = f"RMSD_{ligand}_fitted"
        
        ligand_atom_ids = full_dict[complex_id]["ligands"][ligand]

        if len(ligand_atom_ids) < 3:
            return f"{ligand} is composed of {len(ligand_atom_ids)} atoms, so it does not have self fitted RMSD values."
        else:
            use_superposition = True

        ligand_atom_group = u.atoms[ligand_atom_ids]
        r = rms.RMSD(ligand_atom_group, ligand_atom_group, select = LIGAND_FITTING_CRITERIA, superposition = use_superposition)
        r.run()

        run_df = pd.DataFrame(r.results.rmsd, columns = ["Frame", "Time(ns)"] +  [ligand])
        print(run_df, "\n")
        export_csv_plot(export_name, run_df)

    return "Ligand fitted RMSD successfuly generated. \n"

def compute_atom_distance(u, atoms_ids_listoftuples, max_st):
    """
    Docstring for atom_distances_analysis
    
    :param u: Description
    :param atoms_ids_listoftuples: Description
    """

    group_1_ids = []
    group_2_ids = []
    column_names = []
    for pair in atoms_ids_listoftuples:
        group_1_ids.append(pair[0]-1)
        group_2_ids.append(pair[1]-1)
        column_names.append(f"{pair[0]} / {pair[1]}")
    atom_1_group = u.select_atoms(ids_to_key(group_1_ids))
    atom_2_group = u.select_atoms(ids_to_key(group_2_ids))

    dist_data = []
    for ts in u.trajectory:
        frame = ts.frame
        _, _, atom_distance = distances.dist(atom_1_group, atom_2_group)
        atom_distance_list = [float(atom_distance[i]) for i in range(len(atom_distance))]
        dist_data.append([frame, u.trajectory.time] + atom_distance_list)
        if frame == max_st:
            break
    u.trajectory[0]
    dist_df = pd.DataFrame(dist_data, columns = ["Frame", "Time"] + column_names)

    export_name = "atom_pair_distances"
    export_csv_plot(export_name, dist_df, labels = ["Frames", "Distance (Å)"])
    return "Atom distances successfully generated. \n"



def export_csv_plot(file_export_name, df, labels = None):
    """
    Exports a pandas dataframe as a .csv and a .png (image) files.
    
    :param file_export_name: export name for the .csv and .png files
    :param df: dataframe to be converted
    :param df: list containing optional x and y labels
    """
    colnames = list(df)
    x_axis = colnames[0]
    y_axis = colnames[2:]
    df.to_csv(f"{file_export_name}.csv", index = False)
    df.plot(x = x_axis, y = y_axis)
    plt.legend()
    if labels is not None:
        plt.xlabel(labels[0])
        plt.ylabel(labels[1])
    plt.savefig(f"{file_export_name}.png", dpi = 300)
    plt.close()


def ids_to_key(list):
    out = ""
    for item in list:
        out += f"{item} "
    return f"index {out.strip()}"

def complex_to_ids(full_dict, complex_id):
    """
    Docstring for complex_to_keys
    
    :param full_dict: Description
    :param complex_id: Description
    """
    ids_list= []
    chosen_complex = full_dict[f"{complex_id}"]
    for k in chosen_complex.keys():
        if k == "protein_chain":
            ids_list.append(chosen_complex[k])
        else:
            for kk in chosen_complex[k].keys():
                ids_list.append(chosen_complex[k][kk])
    return ids_list


def get_complex_data(u, full_dict, complex_id):
    """
    From a universe, a dictionary containing all atoms in the system and a specific complex id, \
    outputs the atomids of the complex, the atom_group and the separate keys for all the groups in \
    the complex.
    
    :param full_dict: Description
    :param complex_id: Description
    """
    complex_groups_atomids = complex_to_ids(full_dict, complex_id)

    all_complex_atomids = []
    for group_list in complex_groups_atomids:
        for atom in group_list:
            all_complex_atomids.append(atom)
    complex_atom_group = u.atoms[all_complex_atomids]

    complex_keys_list = []
    for group_atomids in complex_groups_atomids:
        group_keys = ids_to_key(group_atomids)
        complex_keys_list.append(group_keys)

    return complex_atom_group, complex_keys_list



def main():
    RUN_ONCE = ["Distance between atoms"]

    #create universe
    U = mda.Universe(sys.argv[1],sys.argv[2])
       
    #select the mobile domain of the protein specifically
    protein_full = U.select_atoms(PROTEIN_KEY)
    #analyse the radius of gyration
    protein_dict = parse_protein(protein_full)
    #select the ligands specifically
    full_dict = add_ligands_to_complexes(U, protein_dict, LIGANDS_KEYS)

    #compute the RMSD for each complex
    #all functions besides atom distance are run for each complex
    ligand_names = list(full_dict["1"]["ligands"].keys())

    functions_dict = {"RMSD fitted ligands": compute_RMSD_fit_protein ,
                      "RMSD fitted protein": compute_RMSD_fit_ligand ,
                      "Distance between atoms": compute_atom_distance
    }

    for function_name in FUNCTIONS:
        for complex_id in list(full_dict.keys()):

            inputs_dict = {"RMSD fitted ligands": (U, full_dict, complex_id, ligand_names),
                           "RMSD fitted protein": (U, full_dict, complex_id, ligand_names),
                           "Distance between atoms": (U, ATOM_PAIRS, STOP_TIME)
            }
            current_function = functions_dict[function_name]
            current_input = inputs_dict[function_name]
            current_function(*current_input)
            if function_name in RUN_ONCE:
                break


    
    
    


if __name__ == "__main__":
    main()

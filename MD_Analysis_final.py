import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import calc_bonds
import sys
import pandas as pd
from matplotlib import pyplot as plt

#start and stop frames for the analysis functions
#for default values substitute with None
START_FRAME = None
STOP_FRAME = None

# the protein key is used to select residues in all monomers
# the protein is parsed into different monomers automatically
PROTEIN_KEY = "resid 1-50"
PROTEIN_FITTING_CRITERIA = "backbone"

# if more ligands are needed, add them using the reffered notation
# specific atom IDs can selected for the ligands 
LIGANDS_KEYS = ["resname SAM", "resname NAS"]
LIGAND_FITTING_CRITERIA = "all" # "not name H*" to exclude ligand hydrogens

#atom pair indices
#for atom distance analysis
ATOM_PAIRS = [(49, 50), (99, 100), (101, 103)]


#write the functions you want to use
FUNCTIONS = ["Distance between atoms", "Radius of gyration", "RMSD fitted protein", "RMSD fitted ligands"]

"""
Choose functions from function selection pool.
Functions are applied from left to right function options:

RMSD fitted ligands - calculates the RMSD over all frames for each inputed ligand,
fitted for each ligand

RMSD fitted protein - calculates the RMSD over all frames for the protein and 
the ligands, fitted for the inputted protein backbone

Distance between atoms - calculates the distance between 2 atoms over all frames
-> for more that an atom pair, input as a list of tuples: [(A1, A2), (A3, A4)]

Radius of gyration - calculates the radius of gyration (rgyr) over all frames
for each inputed ligand
"""


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
    Adds the ligands of each chain complex in the system to a pre-existing complexes dict
    For instance with 2 complexes and 4 ligands:
    (ligand 1 -> complex 1)
    (ligand 2 -> complex 2) 
    (ligand 3 -> complex 1)
    (ligand 4 -> complex 2)
    
    :param u: Universe
    :param complexes_dict: Complexes dictionary that contains protein chains for each chain complex
    :param ligand_key_list: List contaning ligand keys. Ex: "resname SAM"
    """
    final_dict = complexes_dict.copy()
    all_complex_keys = list(final_dict.keys())
    n_complexes = len(all_complex_keys)

    for ligand_key in ligand_key_list:
        ligand_residues = u.select_atoms(ligand_key).residues
        
        for i, ligand_res in enumerate(ligand_residues):

            ligand_name = ligand_res.resname
            ligand_res_atom_ids = ligand_res.atoms.indices.tolist()
            #distributes the according to the total number of complexes
            #complex names start on 1
            correct_target_id = i % n_complexes
            correct_complex_key = all_complex_keys[correct_target_id]


            if "ligands" not in final_dict[correct_complex_key]:
                final_dict[correct_complex_key]["ligands"] = {}
            
            if ligand_name not in final_dict[correct_complex_key]["ligands"]:
                final_dict[correct_complex_key]["ligands"][ligand_name] = []

            final_dict[correct_complex_key]["ligands"][ligand_name].extend(ligand_res_atom_ids)
        
    return final_dict


def compute_RMSD_fit_protein(u, full_dict, complex_id, ligand_names, protein_fc, min_frame = 0, max_frame = None):
    """
    Computes the RMSD of a inputed complex with fitting on the complex protein_chain, from a universe and a dictionary of all atoms in the system
 
    :param u: Universe
    :param full_dict: Dictionary containing atom ids, separated into complexes. Each complex contains \
    the specific monomer and ligands atoms.
    :param complex_n: Identifier of a complex
    :param protein_fc: The protein fitting criteria used to calculate RMSD values
    :min_frame: Starting frame for the analysis
    :max_frame: Ending frame for the analysis
    """
    first_frame = min_frame
    final_frame = max_frame or len(u.trajectory)
    
    complex_atom_group, complex_keys_list = get_complex_data(u, full_dict, complex_id)

    r = rms.RMSD(complex_atom_group, complex_atom_group, select = protein_fc, groupselections = complex_keys_list)
    r.run(start = first_frame, stop = final_frame)

    run_df = pd.DataFrame(r.results.rmsd, columns = ["Frame", "Time(ns)", "Backbone", protein_fc] +  ligand_names)
    
    export_name = f"RMSD_Protein_complex_{complex_id}_fitted_plus_{'_'.join(ligand_names)}"
    export_csv_plot(export_name, run_df, labels = ["Frames", "Distance (Å)"])

    return "Protein fitted RMSD successfuly generated. \n"


def compute_RMSD_fit_ligand(u, full_dict, complex_id, ligand_names, ligand_fc, min_frame = 0, max_frame = None):
    """
    Computes the RMSD of a inputed complex with fitting on the ligand, from a universe and a dictionary of all atoms in the system
    
    :param u: Universe
    :param full_dict: Dictionary containing atom ids, separated into complexes. Each complex contains \
    the specific monomer and ligands atoms.
    :param complex_id: Identifier of a complex
    :param ligand_names: List containing the ligand names that are to be analysed
    :param ligand_fc: The ligand fitting criteria used to calculate RMSD values
    :min_frame: Starting frame for the analysis
    :max_frame: Ending frame for the analysis
    """

    first_frame = min_frame
    final_frame = max_frame or len(u.trajectory)

    for ligand in ligand_names:
        
        ligand_atom_ids = full_dict[complex_id]["ligands"][ligand]

        if len(ligand_atom_ids) < 3:
            print(f"{ligand} is composed of {len(ligand_atom_ids)} atoms, so it does not have self fitted RMSD values.")
            continue

        ligand_atom_group = u.atoms[ligand_atom_ids]
        r = rms.RMSD(ligand_atom_group, ligand_atom_group, select = ligand_fc)
        r.run(start = first_frame, stop = final_frame)
      
        run_df = pd.DataFrame(r.results.rmsd, columns = ["Frame", "Time(ns)"] +  [ligand])
        
        export_name = f"RMSD_{ligand}_complex_{complex_id}_fitted"
        print(run_df, "\n")
        export_csv_plot(export_name, run_df, labels = ["Frames", "Distance (Å)"])

    return "Ligand fitted RMSD successfuly generated. \n"


def compute_atom_distance(u, atoms_ids_listoftuples, min_frame = 0, max_frame = None):
    """
    Calculates and plots the difference between atom pairs. Accepts the a list of tuples as input for the atom pairs and \
    can calculate the distance for however many atom pairs are included.  

    :param u: Universe
    :param atoms_ids_listoftuples: List of tuples containing the atom pairs. Ex: [(atom_1_from_pair_A, atom_2_from_pair_A), (atom_1_from_pair_B, atom_2_from_pair_B)]
    :min_frame: Starting frame for the analysis
    :max_frame: Ending frame for the analysis
    """

    group_1_ids = []
    group_2_ids = []
    column_names = []
    for pair in atoms_ids_listoftuples:
        group_1_ids.append(pair[0]-1)
        group_2_ids.append(pair[1]-1)
        column_names.append(f"{pair[0]} / {pair[1]}")
    atom_1_group = u.atoms[group_1_ids]
    atom_2_group = u.atoms[group_2_ids]

    first_frame = min_frame
    last_frame = max_frame or len(u.trajectory)

    dist_data = []
    for ts in u.trajectory[first_frame:last_frame]:
        frame = ts.frame
        atom_distance = calc_bonds(atom_1_group.positions, atom_2_group.positions, box=u.dimensions)
        atom_distance_list = [float(atom_distance[i]) for i in range(len(atom_distance))]
        dist_data.append([frame, u.trajectory.time] + atom_distance_list)

    u.trajectory[0]
    dist_df = pd.DataFrame(dist_data, columns = ["Frame", "Time"] + column_names)

    export_name = "atom_pair_distances"
    export_csv_plot(export_name, dist_df, labels = ["Frames", "Distance (Å)"], ind_start= 2)
    return "Atom distances successfully generated. \n"


def analyse_rgyr(u, full_dict, complex_id, ligand_names, min_frame = 0, max_frame = None):
    """
    Calculates the radius of giration for all ligands present in ligand names.
    
    :param u: Universe
    :param full_dict: Dictionary including all atom indexes in the universe
    :param complex_id: Identifier for the complex
    :param ligand_names: List containing the names of all ligands to be analysed
    :min_frame: Starting frame for the analysis
    :max_frame: Ending frame for the analysis
    """
    u.trajectory[0].frame
    first_frame = min_frame
    last_frame = max_frame or len(u.trajectory)

    for ligand in ligand_names:
        
        ligand_atom_ids = full_dict[complex_id]["ligands"][ligand]  
        ligand_atom_list = u.atoms[ligand_atom_ids] 

        rgyr_data = []
        for ts in u.trajectory[first_frame:last_frame]:
            frame = ts.frame
            rgyr_data.append((frame, u.trajectory.time, ligand_atom_list.radius_of_gyration()))
        
        u.trajectory[0]
        rgyr_df = pd.DataFrame(rgyr_data, columns = ["Frame", "Time", "Rgyr"])
        
        export_name = f"Rgyr_{ligand}_complex_{complex_id}"        
        export_csv_plot(export_name, rgyr_df,labels = ["Frames", "Distance (Å)"], ind_start=2)
        return "Radius of gyration analysis sucessfully generated. \n"


def export_csv_plot(file_export_name, df, labels = None, ind_start = 3):
    """
    Exports a pandas dataframe as a .csv and a .png (image) files.
    
    :param file_export_name: Export name for the .csv and .png files
    :param df: Dataframe used for data conversion
    :param df: List containing optional x and y labels
    :param ind_start: Optional input for starting index for y axis values, retrieved \
    from the dataframe colnames 
    """
    colnames = list(df)
    x_axis = colnames[0]
    y_axis = colnames[ind_start:]
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
    Outputs all atoms that belong to the inputed complex_id, from full_dict
    
    :param full_dict: Dictionary including all atom indexes in the universe
    :param complex_id: Identifier for the complex
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
    
    :param u: 
    :param full_dict: Dictionary including all atom indexes in the universe
    :param complex_id: Identifier for the complex
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

    # defines what functions must only run once per command
    RUN_ONCE = ["Distance between atoms"]

    global FUNCTIONS
    FUNCTIONS = list(dict.fromkeys(FUNCTIONS))


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
    

    functions_dict = {"RMSD fitted protein": compute_RMSD_fit_protein ,
                      "RMSD fitted ligands": compute_RMSD_fit_ligand ,
                      "Distance between atoms": compute_atom_distance,
                      "Radius of gyration": analyse_rgyr
    }

    for function_name in FUNCTIONS:
        for complex_id in list(full_dict.keys()):
            ligand_names = list(full_dict[complex_id]["ligands"].keys())

            inputs_dict = {"RMSD fitted ligands": (U, full_dict, complex_id, ligand_names, LIGAND_FITTING_CRITERIA, START_FRAME, STOP_FRAME),
                           "RMSD fitted protein": (U, full_dict, complex_id, ligand_names, PROTEIN_FITTING_CRITERIA, START_FRAME, STOP_FRAME),
                           "Distance between atoms": (U, ATOM_PAIRS, START_FRAME, STOP_FRAME),
                           "Radius of gyration": (U, full_dict, complex_id, ligand_names, START_FRAME, STOP_FRAME)
            }
            current_function = functions_dict[function_name]
            current_input = inputs_dict[function_name]
            current_function(*current_input)
            if function_name in RUN_ONCE:
                break


if __name__ == "__main__":
    main()

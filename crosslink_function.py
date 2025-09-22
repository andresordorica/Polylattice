import os
import mbuild as mb 
from parmed import load_file
import re
from mbuild.lib.atoms import H
import numpy as np
import random
import pubchempy as pcp
from rdkit import Chem
import requests
from collections import Counter


def normalize_crosslinker_input(crosslinker_backbone_target):
    """
    Normalize crosslinker_backbone_target so it's always a list of dictionaries.
    """
    if isinstance(crosslinker_backbone_target, dict):
        return [crosslinker_backbone_target]
    elif isinstance(crosslinker_backbone_target, list):
        if all(isinstance(item, dict) for item in crosslinker_backbone_target):
            return crosslinker_backbone_target
        else:
            raise ValueError("List must contain only dictionaries.")
    else:
        raise TypeError("crosslinker_backbone_target must be a dict or list of dicts.")


def is_dict_of_dicts(d):
    return all(isinstance(v, dict) for v in d.values())


def print_all_pubchem_properties(cid_or_name):
    try:
        # Automatically handle CID or name input
        if isinstance(cid_or_name, int):
            compound = pcp.Compound.from_cid(cid_or_name)
        else:
            compound = pcp.get_compounds(cid_or_name, 'name')[0]

        # Convert to dictionary with all properties
        compound_dict  = compound.to_dict()

        # Print all properties
        #print(f"All properties for: {compound.iupac_name or compound.synonyms[0]}")
        

    except Exception as e:
        print(f"Error retrieving compound: {e}")
    return compound_dict


def extract_absolute_smiles(record):
        if "record" not in record:
            return None

        for entry in record["record"].get("props", []):
            urn = entry.get("urn", {})
            if urn.get("label") == "SMILES":
                value = entry.get("value", {})
                smiles = value.get("sval")
                if smiles:
                    return smiles
        return None
    

def is_match(list_names, bonded_neighbors_copy):
    counter_list = Counter(list_names)
    counter_bonded = Counter(bonded_neighbors_copy)
    
    for elem, count in counter_bonded.items():
        if counter_list[elem] < count:
            return False
    return True
################################### Support functions #######################################################################################################################################
###################################                    ######################################################################################################################################

def fetch_and_plug_CR(name_of_compound):
    '''
    This compound takes a string
    Goes to the pubChem structure database and retrieves the CID (identifier)
    Then it returns the SMILES (isomeric) and the CID of the corresponding molecule
    '''
    result = (pcp.get_compounds(name_of_compound, 'name'))[0]
    digits = result.isomeric_smiles
    cid = result.cid
    print("Base name: {}, SMILES: {}, CID: {}".format(name_of_compound, digits, str(cid)))
    return digits, cid

def get_mol_CR(CID):
    '''
    this function takes the CID of a compounds
    uses RDKIT to determine the mol representation
    can be drawn, plottet etc,
    '''
    '''
    c = pcp.Compound.from_cid(CID)
    print(c)
    smiles=c.isomeric_smiles
    '''    
    smiles = pcp.get_compounds(CID, 'cid', properties=['isomeric_smiles'])[0]    
    record = print_all_pubchem_properties(CID)  # or 'malic acid'  
    smiles = extract_absolute_smiles(record)
    print(smiles)
    mol = Chem.MolFromSmiles(smiles)
    return smiles, mol

def Build_3D_CR(compound_smiles):
    compound = mb.load(compound_smiles,smiles = True, backend = "pybel")
    return compound

def compound_info_CR(parmed_structure):
    atom_info_dict = {}

    for atom in parmed_structure.atoms:
        atom_name = atom.name
        atom_type = atom.type

        # Check if atom_name is already a key in the dictionary
        if atom_name in atom_info_dict:
            # If it is, add the atom_type to the set of atom_types
            atom_info_dict[atom_name].add(atom_type)
        else:
            # If it's not, create a new entry in the dictionary with atom_name as the key
            # Initialize a set with the first atom_type encountered
            atom_info_dict[atom_name] = {atom_type}

    # Convert sets to lists if needed
    for atom_name, atom_types in atom_info_dict.items():
        atom_info_dict[atom_name] = list(atom_types)
        
    return atom_info_dict
        
def Add_port_CR(compound, type = "OH", name = "side",  factor = 1):

    type_anchor = { "OH": [["O", "H"], ["H"], [2], ],
                   "OHt": [["O", "H"], ["H"], [1], ],
                   "NH3": [["N", "H", "H"], ["H"],  [3],], ### Terminal one
                   "NH3t" : [["N", "H"], ["H"],  [2],], ### Terminal one
                   "NH2": [["N", "H", "H"], ["H"],  [3],],
                   "NH": [["N", "C", "H"], ["H"],  [3],],
                   "SH": [["S","H"], ["H"],  [2],],
                   "CH3": [["C", "H", "H", "H"], ["H"], [4]],
                   "CH2": [["C", "C", "H", "H"], ["H"], [3]],
                   "COH": [["C", "C", "O", "H", "H"], ["H"], [4]],
                   "Aryl-alkyl": [["C", "C", "C" ,"H"], ["H"],  [3]],
                   "C=O": [["C", "O", "C" ], ["H"],  [3]],
                   "S=O" :  [["S", "O" ], ["O"],  [3]],
                   "NHc" :  [["N", "H", "H", "C" ], ["H"],  [3]],
                   "CRL_Phenyl" :  [["N", "H", "H", "C" ], ["H"],  [3]],
                   "CRL_Site" :  [["C", "C", "B", "H" ], ["H"],  [3]],
                   "NH2_C": [["N", "H", "H", "C"], ["H"],  [3],],
                   "CO_g": [["C", "O",], ["O"],  [2],],
                   
                   "DGPA_HO": [["B", "H", "C",], ["H"],  [2],],
                   
                   "DGPA_HC": [["C", "H", "C", "O",], ["H"], [3],], 


                    }
    anchor = ((type_anchor[type])[0])[0]
    matchers = ((type_anchor[type])[0])[1:]
    target = ((type_anchor[type])[1])[0]
    nbonds = ((type_anchor[type])[2])[0]

    co =0
    circle = 0
    match = False
    for particle in compound.particles():
        ######################## Case OH
        circle+=1
        if co ==0:
            if particle.name =="{}".format(anchor) and particle.n_direct_bonds ==nbonds:
                list_match = []
                pos_list = []
                part_id = []
                for a in particle.direct_bonds():
                    part_id.append(a)
                    list_match.append(a.name)
                    pos_list.append(a.pos)
                ###########################################################
                from collections import Counter
                # Count the occurrences of each element in both lists
                count1 = Counter(matchers)
                count2 = Counter(list_match)
                var = False
                if type == "Aryl-alkyl":
                    if circle ==5 or circle ==1:
                        var = True
                    else:
                        var = False
                else:
                    var = True
                if var ==True:
                    # Check if all elements in list1 are in list2 with the same or greater counts
                    if all(count1[element] <= count2[element] for element in count1):
                        #print("The matched elements match")
                        position_H = pos_list[list_match.index(target)]
                        import numpy as np
                        # Define the positions of the two points
                        point1 = position_H
                        point2 = particle.pos
                        # Calculate the vector from point1 to point2
                        vector = point1 - point2
                        # Calculate the magnitude (length) of the vector
                        magnitude = np.linalg.norm(vector)
                        # Calculate the unit vector (normalize the vector)
                        unit_vector = vector / magnitude
                        distance = np.linalg.norm(point2 - point1)
                        compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * factor), name)
                        compound.remove(part_id[list_match.index(target)])
                        co+=1
                        match = True
                    else:
                        #print("No match")
                        print("")

                else:
                    print("")
                                 
    return compound

#############################################################################################################################################################################################


def Builder_pipeline_Polymer(name_of_compound, first_port= "OH", second_port = "CH2", add_ports = True, rigid_port = False, fetch_compound = True, MB_compound = None, factor_length = 1):
    print("Polymer loaded as mb compound, no query needed")
    if fetch_compound == True:
        digits, cid = fetch_and_plug_CR(name_of_compound)
        smiles, mol = get_mol_CR(cid)
        record = print_all_pubchem_properties(cid) 
        smiles = extract_absolute_smiles(record)
        compound_mol2 = Build_3D_CR(smiles)
        print("Built the molcule: {} with smiles {}".format(name_of_compound, digits))
        print("Canonical smiles: {} CID:  {}".format(smiles, cid))
    else:
        compound_mol2 = MB_compound
        
    if rigid_port == False: 
        anchor = first_port
        print(anchor)
        print("The port is type: {}".format(anchor))
        if add_ports ==True:
            compound_1_port = Add_port_CR(compound_mol2, type = anchor, name = "up", factor = factor_length)
        else:
            compound_1_port = compound_mol2
        print("Added ports to the molecule")
        base = mb.clone(compound_1_port)

        if add_ports ==True:
            compound_with_2_ports = Add_port_CR(base, type = second_port, name = "down", factor =factor_length )
        else:
            compound_with_2_ports = base
        print("Added 2 ports to the molecule")
     

    return compound_with_2_ports, compound_1_port



def Build_crosslinker(name_of_compound, first_port= "OH", add_ports = True, name_port = "site_1",  fetch_compound = True, MB_compound = None, factor_length = 1):

    if fetch_compound == True:
        digits, cid = fetch_and_plug_CR(name_of_compound)
        smiles, mol = get_mol_CR(cid)
        compound_mol2 = Build_3D_CR(smiles)
        print("Built the molcule: {} with smiles {}".format(name_of_compound, digits))
        print("Canonical smiles: {} CID:  {}".format(smiles, cid))
    else:
        compound_mol2 = MB_compound
            
            
    anchor = first_port
    print(anchor)
    print("The port is type: {}".format(anchor))
    if add_ports ==True:
        compound_1_port = Add_port_CR(compound_mol2, type = anchor, name = name_port, factor = factor_length)
    else:
        compound_1_port = compound_mol2
    print("Added ports to the molecule")


       
        
    return  compound_1_port

class CH2(mb.Compound):
    def __init__(self):
        super(CH2, self).__init__()
        carbon = mb.Particle(pos=[0.0, 0.0, 0.0], name='C')
        hydrogen0 = mb.Particle(pos=[0.109, 0.0, 0.0], name='H')
        hydrogen1 = mb.Particle(pos=[-0.109, 0.0, 0.0], name='H')
        self.add([carbon, hydrogen0, hydrogen1])
        self.add_bond((carbon, hydrogen0))
        self.add_bond((carbon, hydrogen1))
        self.add(mb.Port(anchor=self[0], orientation=[0, 1, 0],
                         separation=0.07), label='up')
        self.add(mb.Port(anchor=self[0], orientation=[0, 0, 1],
                         separation=0.07), label='down')

class OH(mb.Compound):
    def __init__(self):
        super(OH, self).__init__()
        self.add(mb.Particle(name='O', pos=[0.0, 0.0, 0.0]), label='O')
        self.add(mb.Particle(name='H', pos=[0.0, 0.1, 0.0]), label='H')
        self.add_bond((self['O'], self['H']))
        # add the port to the oxygen atom along the [0,-1, 0] direction
        self.add(mb.Port(anchor=self['O'], orientation=[0, -1, 0], separation=0.075), label='down')


class Alcohol(mb.Compound):
    def __init__(self, chain_length=1, up_cap=None):
        super(Alcohol, self).__init__()

        # check if nothing was passed as up_cap
        if up_cap is None:
            # input default behavior here
            up_cap = OH()
        # something was passed as up_cap,
        else:
            pass
        
        # Add to our `self` compound
        self.add(up_cap)
        
        # create the backbone using polymer
        internal_chain = mb.recipes.Polymer(monomers=[CH2()])
        internal_chain.build(n=chain_length, add_hydrogens=False)

        
        # top capping CH2
        mb.force_overlap(move_this=up_cap,
                         from_positions=up_cap['down'],
                         to_positions=internal_chain['up'])

        self.add(internal_chain)

        # bottom cap
        hydrogen = H()
        mb.force_overlap(move_this=hydrogen,
                         from_positions=hydrogen['up'],
                         to_positions=internal_chain['down'])
        self.add(hydrogen)
        
        
def Polymer_anchor_atoms(Polymer_compound = mb.Compound(),
                         Target_atom_positions = dict(),
                         polymer_backbone_targets = list(),
                         bonded_neighbors_bool = True,
                         
                               ):
    list_target_ids = []
    
    c = 0
    for key, value in Target_atom_positions.items():
        for child in Polymer_compound.children:
            if key == child.name:
                for dictionary_targets in polymer_backbone_targets:
                    anchor_atom = dictionary_targets['anchor_atom']
                    anchor_atom_n_bonds = dictionary_targets['anchor_atom_n_bonds']
                    bonded_neighbors = dictionary_targets['bonded_atoms']
                    bonded_neighbors_copy = bonded_neighbors.copy()
                    if anchor_atom in bonded_neighbors_copy:
                        bonded_neighbors_copy.remove(anchor_atom)
                        
                    for particle in child.particles():
                        if particle.name ==anchor_atom:
                            list_names = []
                            list_positions = []
                            list_id = []
                                                                
                            for p in particle.direct_bonds():
                                list_names.append(p.name)
                                list_positions.append(p.pos)
                                list_id.append(p)
                                
                            if len(list_names) == anchor_atom_n_bonds:
                                # Convert to counters
                                if bonded_neighbors_bool:
                                    #print(f"anchor atom {anchor_atom}")
                                    #print(f"List of neighbos: {list_names}, list dict {bonded_neighbors_copy}")
                                    if is_match(list_names, bonded_neighbors_copy):
                                        #print("match")
                                        Target_atom_positions[key].append((c, particle.pos))
                                        list_target_ids.append(id(particle))
                                        c += 1
                                    
                                else:
                                    Target_atom_positions[key].append((c, particle.pos))
                                    list_target_ids.append(id(particle))
                                    c += 1
                                
                            
          
            
    return Target_atom_positions, list_target_ids
                               

def select_continuous_groups_random(data, key_name, N_number_of_groups, len_of_groups, number_of_polymers, force_even_start=False, tuple_case = True):
    if key_name not in data:
        raise ValueError(f"Key '{key_name}' not found in data.")

    tuples_list = data[key_name]  # Retrieve the list associated with the key
    total_tuples = len(tuples_list)

    
    if tuple_case == True:
        # Ensure we have enough tuples to satisfy the request
        if number_of_polymers <=2:
            required_tuples = (N_number_of_groups * len_of_groups)//2
        else:
            required_tuples = N_number_of_groups * len_of_groups
            
        if required_tuples > total_tuples:
            raise ValueError(f"Not enough tuples in '{key_name}'. "
                            f"Requested {required_tuples}, but only {total_tuples} available.")

        # Generate valid starting indices
        if force_even_start:
            available_indices = [i for i in range(total_tuples - len_of_groups + 1) if i % 2 == 0]
        else:
            available_indices = list(range(total_tuples - len_of_groups + 1))

        random.shuffle(available_indices)  # Shuffle to introduce randomness

        selected_dict = {}
        used_indices = set()  # Keep track of used indices to avoid repetition

        for i in available_indices:
            if len(selected_dict) >= N_number_of_groups:
                break  # Stop when we have enough entries

            group = tuples_list[i:i + len_of_groups]
            group_indices = {t[0] for t in group}  # Extract indices from the tuples
            #print(group)
            #print(group_indices)
            # Ensure the group is continuous and indices have not been used
            if all(group[j][0] + 1 == group[j + 1][0] for j in range(len_of_groups - 1)) and not used_indices.intersection(group_indices):
                for index, position in group:
                    selected_dict[index] = position  # Assign to dictionary
                used_indices.update(group_indices)  # Mark these indices as used
                
    else:
        # Ensure we have enough tuples to satisfy the request
        if number_of_polymers <=2:
            required_tuples = (N_number_of_groups * len_of_groups)//2
        else:
            required_tuples = N_number_of_groups * len_of_groups
            
        if required_tuples > total_tuples:
            raise ValueError(f"Not enough tuples in '{key_name}'. "
                            f"Requested {required_tuples}, but only {total_tuples} available.")

        # Generate valid starting indices
        #if force_even_start:
            #available_indices = [i for i in range(total_tuples - len_of_groups + 1) if i % 2 == 0]
        #else:
        available_indices = list(range(total_tuples - len_of_groups + 1))

        random.shuffle(available_indices)  # Shuffle to introduce randomness

        selected_dict = {}
        used_indices = set()  # Keep track of used indices to avoid repetition

        for i in available_indices:
            if len(selected_dict) >= N_number_of_groups:
                break  # Stop when we have enough entries

            group = tuples_list[i]
            group_indices = {group[0]}  # Extract indices from the tuples
            #print(group)
            #print(group_indices)
            # Ensure the group is continuous and indices have not been used
            if not used_indices.intersection(group_indices):
                selected_dict[group[0]] = group[1]  # Assign to dictionary
                used_indices.update(group_indices)  # Mark these indices as used
        
    #print("Selected dictionary")
    #print(selected_dict)
    return selected_dict
 

def check_position_in_dict(position_array, position_dict, tolerance=1e-2):
    """
    Checks if a given position array matches any of the values in the dictionary within a given tolerance.

    Parameters:
    - position_array: The array to check.
    - position_dict: Dictionary where values are position arrays.
    - tolerance: The acceptable difference for a match (default: 1e-2).

    Returns:
    - True if a match is found within the tolerance, False otherwise.
    """
    for stored_array in position_dict.values():
        if np.allclose(position_array, stored_array, atol=tolerance):
            return True
        
    return False


def extract_matching_tuples(data, key_name,):
    """
    Extracts tuples from a given dictionary where the first element of the tuples matches values in base_list.

    Parameters:
    - data: Dictionary containing lists of tuples.
    - key_name: The key in the dictionary to extract from.
    - base_list: List of values to match with the first element of the tuples.

    Returns:
    - A dictionary where matched first elements are keys and corresponding position arrays are values.
    """
    
    tuples_list = data[key_name]  # Retrieve the list of tuples
    result_dict = {index: position for index, position in tuples_list}

    return result_dict


def process_dicts(previous_dict, groups_pol_1, end_cap=False, tuple_case = True):
    """
    Removes keys from groups_pol_1 whose first two elements in values match those in previous_dict.
    Then selects all consecutive key pairs from the filtered dictionary without repeating elements.

    Parameters:
    - previous_dict (dict): Dictionary with keys mapping to numpy arrays.
    - groups_pol_1 (dict): Dictionary with keys mapping to numpy arrays.
    - end_cap (bool): If True, returns only pairs that exist in previous_dict and the selected pairs dict
                      contains only elements in common between both dictionaries.

    Returns:
    - removed_keys (int): Number of keys removed.
    - selected_pairs (list): List of tuples containing all consecutive key pairs.
    - selected_pairs_dict (dict): Dictionary with selected key pairs and their corresponding values.
    """
    # Extract first two elements of each value in previous_dict
    previous_pairs = {tuple(v[:2]) for v in previous_dict.values()}

    # Filter groups_pol_1 to remove matching keys
    filtered_groups_pol_1 = {k: v for k, v in groups_pol_1.items() if tuple(v[:2]) not in previous_pairs}

    # Count removed keys
    removed_keys = len(groups_pol_1) - len(filtered_groups_pol_1)

    if end_cap:
        # Filter previous_dict to only keep keys that exist in groups_pol_1
        filtered_previous_dict = {k: v for k, v in previous_dict.items() if k in groups_pol_1}

        # Sort keys of filtered previous_dict
        sorted_keys = sorted(filtered_previous_dict.keys())

        if tuple_case == True:
            # Find consecutive pairs in previous_dict
            consecutive_pairs = [
                (sorted_keys[i], sorted_keys[i + 1])
                for i in range(len(sorted_keys) - 1)
                if sorted_keys[i] + 1 == sorted_keys[i + 1]
            ]
            selected_pairs_dict = {k: filtered_groups_pol_1[k] for pair in consecutive_pairs for k in pair}
            
        else:
            consecutive_pairs = [
                (sorted_keys[i])
                for i in range(len(sorted_keys))
    
            ]
            #print(filtered_groups_pol_1)
            #print(consecutive_pairs)
            selected_pairs_dict = {}
            for element in consecutive_pairs:
                selected_pairs_dict[element] = filtered_groups_pol_1[element]

            #print(selected_pairs_dict)
            

        # Create a dictionary containing only selected key pairs and their values from previous_dict
        #selected_pairs_dict = {k: previous_dict[k] for pair in consecutive_pairs for k in pair}
    
    else:
        # Sort keys of filtered dictionary
        sorted_keys = sorted(filtered_groups_pol_1.keys())

        if tuple_case == True:
            # Find consecutive pairs in filtered_groups_pol_1
            consecutive_pairs = [
                (sorted_keys[i], sorted_keys[i + 1])
                for i in range(len(sorted_keys) - 1)
                if sorted_keys[i] + 1 == sorted_keys[i + 1]
            ]
            
            selected_pairs_dict = {k: filtered_groups_pol_1[k] for pair in consecutive_pairs for k in pair}
            
        else:
            consecutive_pairs = [
                (sorted_keys[i])
                for i in range(len(sorted_keys))
    
            ]
        
            #print(filtered_groups_pol_1)
            #print(consecutive_pairs)
            selected_pairs_dict = {}
            for element in consecutive_pairs:
                selected_pairs_dict[element] = filtered_groups_pol_1[element]

            #print(selected_pairs_dict)
    
    #print("Selected pairs dictionary")
    #print(selected_pairs_dict)
    return removed_keys, consecutive_pairs, selected_pairs_dict

def filter_selected_pairs(selected_pairs_dict, attached_ports):
    """
    Removes keys from selected_pairs_dict if their values match (in x, y) with values in attached_ports.

    :param selected_pairs_dict: Dictionary with keys mapping to lists (arrays of at least 2 elements).
    :param attached_ports: Dictionary with keys mapping to NumPy arrays of at least 2 elements.
    :return: Filtered dictionary without matching (x, y) values.
    """
    # Extract (x, y) pairs from attached_ports values
    attached_xy = {tuple(val[:2]) for val in attached_ports.values()}

    # Filter out keys whose first two values (x, y) exist in attached_ports values
    product = {k: v for k, v in selected_pairs_dict.items() if tuple(v[:2]) not in attached_xy}
    #print("Final selected pairs")
    #print(product)
    return product

def filter_non_overlapping_continuous_pairs(data, tuple_case = True):
    """
    Filters out elements that have a continuous key-value pair in the dictionary.
    Ensures that keys are only used once in a pair.

    :param data: Dictionary with integer keys and array values.
    :return: Filtered dictionary with only non-overlapping continuous pairs.
    """
    sorted_keys = sorted(data.keys())
    filtered_keys = set()
    result = {}
    
    if tuple_case == True:
        for i in range(len(sorted_keys) - 1):
            if sorted_keys[i] + 1 == sorted_keys[i + 1] and sorted_keys[i] not in filtered_keys:
                result[sorted_keys[i]] = data[sorted_keys[i]]
                result[sorted_keys[i + 1]] = data[sorted_keys[i + 1]]
                filtered_keys.add(sorted_keys[i + 1])  # Mark this key as used to avoid overlapping
                
    else:
       result = data

    #print("Filtered non-overlapping")
    #print(result)
    return result

def Polymer_Child_Remove_atoms(Polymer_compound = mb.Compound(),
                               child_name = "POL-0",
                               anchor_atom = "O",
                               anchor_atom_n_bonds = 2,
                               target_atom = "H",
                               atoms_to_prune = 1,
                               port_name = "Port",
                               bond_length_factor = 2,
                               dictionary_removal = dict(),
                               dictionary_port_index = dict(),
                               restricted_insertion = False,
                               ports_dict = dict(),
                               
                               ):
    
    list_attached_ports = []
    
    if not restricted_insertion:
        c = 0
        ll_1 =[]
        previous_port = {}
        removed_particle_dictionary =  dict()
        for child in Polymer_compound.children:
            if child.name  == child_name:
                list_alpha = []
                for particle in child.particles():
                    if c < atoms_to_prune:
                        if particle.name ==anchor_atom:
                            list_names = []
                            list_positions = []
                            list_id = []
                            boolean_removal = check_position_in_dict(particle.pos, dictionary_port_index)
                            if boolean_removal == True:      
                                for p in particle.direct_bonds():
                                    list_names.append(p.name)
                                    list_positions.append(p.pos)
                                    list_id.append(p)
                                    
                                if (len(list_names) == anchor_atom_n_bonds) and (target_atom in list_names)   : # Modified line
                                    #print(f"Found target {target_atom} atom")
                                    position_H = list_positions[list_names.index(target_atom)]
                                    # Define the positions of the two points
                                    point1 = position_H
                                    point2 = particle.pos
                                    # Calculate the vector from point1 to point2
                                    vector = point1 - point2
                                    # Calculate the magnitude (length) of the vector
                                    magnitude = np.linalg.norm(vector)
                                    # Calculate the unit vector (normalize the vector)
                                    unit_vector = vector / magnitude
                                    distance = np.linalg.norm(point2 - point1)
                                    ###########
                                    name_ = port_name +"_"+ str(c)
                                    Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                    Polymer_compound.remove(list_id[list_names.index(target_atom)])
                                    list_alpha.append((c,point2,point1))
                                    previous_port[c]= particle.pos
                                    ll_1.append(name_)
                                    #print("Removed particle")
                                    c+=1
                                    
                                    
       
                                
        list_attached_ports.append(ll_1)                                                    
            
   
    else:
        ll_1 = []
        previous_port = {}
        removed_particle_dictionary =  dict()
        list_keys = []
        for child in Polymer_compound.children:
            if child.name  == child_name:
                list_alpha = []
                for particle in child.particles():
                   
                    if particle.name ==anchor_atom:
                        list_names = []
                        list_positions = []
                        list_id = []
                        #boolean_removal = check_position_in_dict(particle.pos, dictionary_port_index)
                        #if boolean_removal == True:      
                        for p in particle.direct_bonds():
                            list_names.append(p.name)
                            list_positions.append(p.pos)
                            list_id.append(p)
                            #print(f"target atom {target_atom}")
                        if (len(list_names) == anchor_atom_n_bonds) and (target_atom in list_names)   : # Modified line
                            #print(list_names)
                            #print(f"Found target {target_atom} atom")
                            position_H = list_positions[list_names.index(target_atom)]
                            # Define the positions of the two points
                            point1 = position_H
                            point2 = particle.pos
                            # Calculate the vector from point1 to point2
                            vector = point1 - point2
                            # Calculate the magnitude (length) of the vector
                            magnitude = np.linalg.norm(vector)
                            # Calculate the unit vector (normalize the vector)
                            unit_vector = vector / magnitude
                            distance = np.linalg.norm(point2 - point1)
                            ###########
                            #print(ports_dict)
                            #print("Restricted insertion")
                            for key, val in ports_dict.items():
                                array2 = ports_dict[key]
                                if np.allclose(particle.pos[:2], array2[:2], atol=1e-2):
                                    #print("Matched")
                                    #print( ports_dict[key])
                                    name_ = port_name +"_"+ str(key)
                                    Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                    Polymer_compound.remove(list_id[list_names.index(target_atom)])
                                    #list_alpha.append((c,point2,point1))
                                    #previous_port[key]= particle.pos
                                    #print("Removed particle")
                                    ll_1.append(name_)
                                    list_keys.append(key)
                                    
        list_attached_ports.append(ll_1)                        
        min_val = min(list_keys)
        max_val = max(list_keys)
        print("List of keys")
        print(list_keys)
        if min_val == 0:  # Case where the list starts at 0
            c = max_val + 1
            atoms_to_prune =   atoms_to_prune 
        else:  # Case where the list starts from a non-zero value
            c = 0
            atoms_to_prune =   atoms_to_prune - len(list_keys) 
                                    
        ll_1 = []
        for child in Polymer_compound.children:
            if child.name  == child_name:
                list_alpha = []
                for particle in child.particles():
                    if c < atoms_to_prune:
                        if particle.name ==anchor_atom:
                            list_names = []
                            list_positions = []
                            list_id = []
                            boolean_removal = check_position_in_dict(particle.pos, dictionary_port_index)
                            if boolean_removal == True:      
                                for p in particle.direct_bonds():
                                    list_names.append(p.name)
                                    list_positions.append(p.pos)
                                    list_id.append(p)
                                    
                                if (len(list_names) == anchor_atom_n_bonds) and (target_atom in list_names) : # Modified line
                                    #print(f"Found target {target_atom} atom")
                                    position_H = list_positions[list_names.index(target_atom)]
                                    # Define the positions of the two points
                                    point1 = position_H
                                    point2 = particle.pos
                                    # Calculate the vector from point1 to point2
                                    vector = point1 - point2
                                    # Calculate the magnitude (length) of the vector
                                    magnitude = np.linalg.norm(vector)
                                    # Calculate the unit vector (normalize the vector)
                                    unit_vector = vector / magnitude
                                    distance = np.linalg.norm(point2 - point1)
                                    ###########
                                    name_ = port_name +"_"+ str(c)
                                    Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                    Polymer_compound.remove(list_id[list_names.index(target_atom)])
                                    #list_alpha.append((c,point2,point1))
                                    previous_port[c]= particle.pos
                                    ll_1.append(name_)
                                    #print("Removed particle")
                                    c+=1
                                                            
               
        list_attached_ports.append(ll_1)
        removed_particle_dictionary[child_name] = list_alpha
                

                    
    return Polymer_compound, removed_particle_dictionary, previous_port, list_attached_ports
    


def Find_Shift_value (
                        Monomer = mb.Compound(),
                        targets_dictionary_list = list(),
                        dimension = 0,
                        ):
    
    dictionary_ = {}
    for dict_target in targets_dictionary_list:
        anchor_atom = dict_target['anchor_atom']
        anchor_atom_n_bonds = dict_target['anchor_atom_n_bonds']
        target_atom = dict_target['target_atom']
        bonded_neighbors = dict_target['bonded_atoms']
        bonded_neighbors_copy = bonded_neighbors.copy()
        ####
        dictionary_[anchor_atom] = []
        if anchor_atom in bonded_neighbors_copy:
            bonded_neighbors_copy.remove(anchor_atom)
                        
        for particle in Monomer.particles():
            if particle.name ==anchor_atom:
                list_names = []
                list_positions = []
                list_id = []
                for p in particle.direct_bonds():
                    list_names.append(p.name)
                    list_positions.append(p.pos)
                    list_id.append(p)
                    
                if (len(list_names) == anchor_atom_n_bonds) :
                    if is_match(list_names, bonded_neighbors_copy):
                        #print("Current anchor particle in monomer")
                        #print(particle.name)
                       # print("Current bonded atoms in monomer")
                       # print(list_names)
                        #print("This are the ones the user defined")
                        #print(bonded_neighbors_copy)
                        dictionary_[anchor_atom].append(particle.pos)
    

    def absolute_column_difference(position_dict):
        # Flatten and stack arrays
        all_positions = [pos for sublist in position_dict.values() for pos in sublist]
        stacked = np.vstack(all_positions)
        
        if stacked.shape[0] != 2:
            raise ValueError("Function only supports exactly 2 vectors for comparison.")
        
        # Column-wise absolute difference
        abs_diff = np.abs(stacked[0] - stacked[1])
        return abs_diff
    
    diff_array = absolute_column_difference(dictionary_)
    
    print("This is the shift values dictionary")
    print(dictionary_)
    print("This is the shift values")
    print(diff_array)
    value = diff_array[dimension]
    print(f"The value for dimension {dimension} is : {value}")

    return value   
                        
                        

def Polymer_Child_Remove_atoms_NESTED(Polymer_compound = mb.Compound(),
                               child_name = "POL-0",
                                atoms_to_prune = 1,
                               targets_dictionary_list = list(),
                               port_name = "Port",
                               bond_length_factor = 2,
                               dictionary_removal = dict(),
                               dictionary_port_index = dict(),
                               restricted_insertion = False,
                               ports_dict = dict(),
                               threshold = 1e-2,
                               
                               ):
    
    
    target_list_purged = [targets_dictionary_list[0]]
   
    list_attached_ports = []
    
    if not restricted_insertion:
        c = 0
        ll_1 =[]
        print(f"Atoms to prune: {atoms_to_prune}")
        previous_port = {}
        removed_particle_dictionary =  dict()
        for IDX, dict_target in enumerate(target_list_purged):
            anchor_atom = dict_target['anchor_atom']
            anchor_atom_n_bonds = dict_target['anchor_atom_n_bonds']
            target_atom = dict_target['target_atom']
            for child in Polymer_compound.children:
                if child.name  == child_name:    
                    list_alpha = []
                    for particle in child.particles():
                        if c < atoms_to_prune:
                            if particle.name ==anchor_atom:
                                list_names = []
                                list_positions = []
                                list_id = []
                                #print(dictionary_port_index)
                                boolean_removal = check_position_in_dict(particle.pos, dictionary_port_index, tolerance= threshold)
                                if boolean_removal == True:
                                    print("Found the site one O")      
                                    for p in particle.direct_bonds():
                                        list_names.append(p.name)
                                        list_positions.append(p.pos)
                                        list_id.append(p)
                                        
                                    if (len(list_names) == anchor_atom_n_bonds) and (target_atom in list_names)   : # Modified line
                                        #print(f"Found target {target_atom} atom")
                                        position_H = list_positions[list_names.index(target_atom)]
                                        # Define the positions of the two points
                                        point1 = position_H
                                        point2 = particle.pos
                                        # Calculate the vector from point1 to point2
                                        vector = point1 - point2
                                        # Calculate the magnitude (length) of the vector
                                        magnitude = np.linalg.norm(vector)
                                        # Calculate the unit vector (normalize the vector)
                                        unit_vector = vector / magnitude
                                        distance = np.linalg.norm(point2 - point1)
                                        ###########

                                        name_ = port_name +"_"+ str(c)
                                        print(f"added Port: {name_}")
                                        Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                        Polymer_compound.remove(list_id[list_names.index(target_atom)])
                                        list_alpha.append((c,point2,point1))
                                        previous_port[c]= particle.pos
                                        ll_1.append(name_)
                                        print(f"Removed particle {c}")
                                        c+=1
                                    
                                    
       
                                
        list_attached_ports.append(ll_1)                                                    
            
   
    else:
        ll_1 = []
        previous_port = {}
        removed_particle_dictionary =  dict()
        list_keys = []
        for child in Polymer_compound.children:
            if child.name  == child_name:
                for IDX, dict_target in enumerate(targets_dictionary_list):
                    anchor_atom = dict_target['anchor_atom']
                    anchor_atom_n_bonds = dict_target['anchor_atom_n_bonds']
                    target_atom = dict_target['target_atom']
                    list_alpha = []
                    for particle in child.particles():
                    
                        if particle.name ==anchor_atom:
                            list_names = []
                            list_positions = []
                            list_id = []
                            #boolean_removal = check_position_in_dict(particle.pos, dictionary_port_index)
                            #if boolean_removal == True:      
                            for p in particle.direct_bonds():
                                list_names.append(p.name)
                                list_positions.append(p.pos)
                                list_id.append(p)
                                #print(f"target atom {target_atom}")
                            if (len(list_names) == anchor_atom_n_bonds) and (target_atom in list_names)   : # Modified line
                                #print(list_names)
                                #print(f"Found target {target_atom} atom")
                                position_H = list_positions[list_names.index(target_atom)]
                                # Define the positions of the two points
                                point1 = position_H
                                point2 = particle.pos
                                # Calculate the vector from point1 to point2
                                vector = point1 - point2
                                # Calculate the magnitude (length) of the vector
                                magnitude = np.linalg.norm(vector)
                                # Calculate the unit vector (normalize the vector)
                                unit_vector = vector / magnitude
                                distance = np.linalg.norm(point2 - point1)
                                ###########
                                #print(ports_dict)
                                #print("Restricted insertion")
                                for key, val in ports_dict.items():
                                    array2 = ports_dict[key]
                                    val_arr = array2[0]
                                    val_particle = particle.pos[0]
                                   
                                      
                                    if np.allclose(val_particle, val_arr, atol=0.001):
                                        #print("Matched")
                                        #print("previous chain port pos:")
                                        #print(val_arr)
                                        #print(f"current particle name {particle.name} and pos")
                                        ##print(val_particle)
                                        #print( ports_dict[key])
                                        #print(f"Added the port on the chain on atom {particle.name}")
                                        name_ = port_name +"_"+ str(key)
                                        print(f"Added Port {name_}")
                                        Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                        Polymer_compound.remove(list_id[list_names.index(target_atom)])
                                        #list_alpha.append((c,point2,point1))
                                        #previous_port[key]= particle.pos
                                        #print("Removed particle")
                                        ll_1.append(name_)
                                        list_keys.append(key)
                                    
        list_attached_ports.append(ll_1)                        
        min_val = min(list_keys)
        max_val = max(list_keys)
        print("List of keys")
        print(list_keys)
        print(f"Counter at {len(list_keys)}")
        print(f"atoms to prune {atoms_to_prune}")
        print(list_keys)
        if min_val == 0:  # Case where the list starts at 0
            c = max_val + 1
            atoms_to_prune =   atoms_to_prune 
        else:  # Case where the list starts from a non-zero value
            c = 0
            atoms_to_prune =   atoms_to_prune - len(list_keys) 


       
                               
        ll_1 = []
        for child in Polymer_compound.children:
            if child.name  == child_name:
                for IDX, dict_target in enumerate(target_list_purged):
                    anchor_atom = dict_target['anchor_atom']
                    anchor_atom_n_bonds = dict_target['anchor_atom_n_bonds']
                    target_atom = dict_target['target_atom']
                    list_alpha = []
                    for particle in child.particles():
                        
                        if c < atoms_to_prune:
                            if particle.name ==anchor_atom:
                                list_names = []
                                list_positions = []
                                list_id = []
                                #boolean_removal = check_position_in_dict(particle.pos, dictionary_port_index, tolerance= threshold)
                                #if boolean_removal == True:      
                                for p in particle.direct_bonds():
                                    list_names.append(p.name)
                                    list_positions.append(p.pos)
                                    list_id.append(p)
                                
                            
                                if (len(list_names) == anchor_atom_n_bonds) and (target_atom in list_names) : # Modified line
                                    position_H = list_positions[list_names.index(target_atom)]
                                    # Define the positions of the two points
                                    point1 = position_H
                                    point2 = particle.pos
                                    # Calculate the vector from point1 to point2
                                    vector = point1 - point2
                                    # Calculate the magnitude (length) of the vector
                                    magnitude = np.linalg.norm(vector)
                                    # Calculate the unit vector (normalize the vector)
                                    unit_vector = vector / magnitude
                                    distance = np.linalg.norm(point2 - point1)
                                    ###########
                                    name_ = port_name +"_"+ str(c)
                                    print(f"Added Port {name_}")
                                    Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                    Polymer_compound.remove(list_id[list_names.index(target_atom)])
                                    #list_alpha.append((c,point2,point1))
                                    previous_port[c]= particle.pos
                                    ll_1.append(name_)
                                    #print("Removed particle")
                                    c+=1
                                                            
               
        list_attached_ports.append(ll_1)
        removed_particle_dictionary[child_name] = list_alpha
                
    print(f"Removed len {c} particles")
                    
    return Polymer_compound, removed_particle_dictionary, previous_port, list_attached_ports
    
      

                        

def add_unique_random_numbers(keys_list, N, repeat_units):
    """Add N unique random numbers (0 to repeat_units) not in keys_list."""
    available_numbers = set(range(repeat_units)) - set(keys_list)
    if len(available_numbers) < N:
        raise ValueError("Not enough unique numbers available to add.")
    
    new_numbers = random.sample(available_numbers, N)
    keys_list.extend(new_numbers)
    return keys_list, new_numbers


def Bond_Order_Tuple(Polymer_residue_name, number_of_polymers, number_repeat_units, crosslinkers_per_pair_of_chains, crosslink_sites_per_molecule, total_crosslinkers_Degree_Crosslinking,Crosslinker_port_name   ):


    Ports_list_polymer = []
    for n in range(number_of_polymers):
        each_chain_ports = []
        if n ==0 or n == (number_of_polymers - 1):
            Mid_chain_true = 1
        else:
            Mid_chain_true = 2
        for p in  range(int(crosslinkers_per_pair_of_chains * crosslink_sites_per_molecule* Mid_chain_true)):   # 3 crosslinkk between chains * 2 hydroxyls per chain = 6
            name_port_it = f"{Polymer_residue_name}-{n}_{p}"
            each_chain_ports.append(name_port_it)
            
            
        Ports_list_polymer.append(each_chain_ports)
        
        
    def generate_even_list(N):
        return [i for i in range(0, 2 * N, 2)]

    tuple_bond_order_ports = []
    for Chain_ in range(len(Ports_list_polymer)-1):
        #print(f"Polymer chain: {Chain_}")
        if Chain_ == 0:
            PLP = Ports_list_polymer     
            len_ports_chain = len(PLP[Chain_])
            sub_slice = (PLP[Chain_+1])[:(len_ports_chain)]
            to_remove = []
            for n in generate_even_list(crosslinkers_per_pair_of_chains ):
                tuple_ =  ( (PLP[Chain_])[n],
                        (PLP[Chain_])[n+1],
                        sub_slice[n],
                        sub_slice[n+1]
                        )
                
                to_remove.extend(tuple_)
                tuple_bond_order_ports.append(tuple_)
                
            # Remove stored elements from the list
            for element in to_remove:
                if element in PLP[Chain_]:  # Check before removing
                    PLP[Chain_].remove(element)
                elif element in PLP[Chain_ + 1]:
                    PLP[Chain_ + 1].remove(element)
                    
                    
        elif Chain_ == (len(Ports_list_polymer)-2) : 
            #print("Last iteration")
            len_ports_chain = len(PLP[Chain_])
            sub_slice = (PLP[Chain_+1])[:(len_ports_chain)]
            to_remove = []
            for n in generate_even_list(crosslinkers_per_pair_of_chains ):
                tuple_ =  ( (PLP[Chain_])[n],
                        (PLP[Chain_])[n+1],
                        sub_slice[n],
                        sub_slice[n+1]
                        )
                
                to_remove.extend(tuple_)
                tuple_bond_order_ports.append(tuple_)
                
            # Remove stored elements from the list
            for element in to_remove:
                if element in PLP[Chain_]:  # Check before removing
                    PLP[Chain_].remove(element)
                elif element in PLP[Chain_ + 1]:
                    PLP[Chain_ + 1].remove(element)
                    
        else:
            len_ports_chain = len(PLP[Chain_])
            sub_slice = (PLP[Chain_+1])[:(len_ports_chain)]
            to_remove = []
            for n in generate_even_list(crosslinkers_per_pair_of_chains ):
                tuple_ =  ( (PLP[Chain_])[n],
                        (PLP[Chain_])[n+1],
                        sub_slice[n],
                        sub_slice[n+1]
                        ) 
                
                to_remove.extend(tuple_)
                tuple_bond_order_ports.append(tuple_)
                
            # Remove stored elements from the list
            for element in to_remove:
                if element in PLP[Chain_]:  # Check before removing
                    PLP[Chain_].remove(element)
                elif element in PLP[Chain_ + 1]:
                    PLP[Chain_ + 1].remove(element)
                    
                
        
    Ports_list_crosslinker_tmp = []
    for n in range(total_crosslinkers_Degree_Crosslinking):
        each_crosslinker_ports = []
        for p in  range(int(crosslink_sites_per_molecule * 2)):   # 2 crosslinking sites per side of crosslinker * 2 
            name_port_it = f"{Crosslinker_port_name}-{n}_{p}"
            each_crosslinker_ports.append(name_port_it)
            
            
        Ports_list_crosslinker_tmp.append(each_crosslinker_ports)
        
    tuple_bond_order_sites = [tuple(sublist) for sublist in Ports_list_crosslinker_tmp]
    
    if len(tuple_bond_order_sites) != len(tuple_bond_order_ports):
        raise ValueError("The number of sites tuples does not match the number of ports tuples")

        

    return tuple_bond_order_ports, tuple_bond_order_sites


def Bond_Order_Non_Tuple(Polymer_residue_name,
                         number_of_polymers,
                         number_repeat_units, crosslinkers_per_pair_of_chains, crosslink_sites_per_molecule, total_crosslinkers_Degree_Crosslinking,Crosslinker_port_name   ):


    Ports_list_polymer = []
    for n in range(number_of_polymers):
        each_chain_ports = []
        if n ==0 or n == (number_of_polymers - 1):
            Mid_chain_true = 1
        else:
            Mid_chain_true = 2
            
        for p in  range(int(crosslinkers_per_pair_of_chains * crosslink_sites_per_molecule* Mid_chain_true)):   # 3 crosslinkk between chains * 2 hydroxyls per chain = 6
            name_port_it = f"{Polymer_residue_name}-{n}_{p}"
            print(name_port_it)
            each_chain_ports.append(name_port_it)
            
            
        Ports_list_polymer.append(each_chain_ports)
        
    print(Ports_list_polymer)
    def generate_even_list(N):
        return [i for i in range(0, N, 1)]

    tuple_bond_order_ports = []
    for Chain_ in range(len(Ports_list_polymer)-1):
        #print(f"Polymer chain: {Chain_}")
        if Chain_ == 0:
            PLP = Ports_list_polymer     
            len_ports_chain = len(PLP[Chain_])
            sub_slice = (PLP[Chain_+1])[:(len_ports_chain)]
            to_remove = []
            for n in generate_even_list(crosslinkers_per_pair_of_chains ):
                tuple_ =  ( (PLP[Chain_])[n],
                        #(PLP[Chain_])[n+1],
                        sub_slice[n],
                        #sub_slice[n+1]
                        )
                
                to_remove.extend(tuple_)
                tuple_bond_order_ports.append(tuple_)
                
            # Remove stored elements from the list
            for element in to_remove:
                if element in PLP[Chain_]:  # Check before removing
                    PLP[Chain_].remove(element)
                elif element in PLP[Chain_ + 1]:
                    PLP[Chain_ + 1].remove(element)
                    
                    
        elif Chain_ == (len(Ports_list_polymer)-2) : 
            #print("Last iteration")
            len_ports_chain = len(PLP[Chain_])
            sub_slice = (PLP[Chain_+1])[:(len_ports_chain)]
            to_remove = []
            for n in generate_even_list(crosslinkers_per_pair_of_chains ):
                tuple_ =  ( (PLP[Chain_])[n],
                        #(PLP[Chain_])[n+1],
                        sub_slice[n],
                        #sub_slice[n+1]
                        )
                
                to_remove.extend(tuple_)
                tuple_bond_order_ports.append(tuple_)
                
            # Remove stored elements from the list
            for element in to_remove:
                if element in PLP[Chain_]:  # Check before removing
                    PLP[Chain_].remove(element)
                elif element in PLP[Chain_ + 1]:
                    PLP[Chain_ + 1].remove(element)
                    
        else:
            len_ports_chain = len(PLP[Chain_])
            sub_slice = (PLP[Chain_+1])[:(len_ports_chain)]
            to_remove = []
            for n in generate_even_list(crosslinkers_per_pair_of_chains ):
                tuple_ =  ( (PLP[Chain_])[n],
                       # (PLP[Chain_])[n+1],
                        sub_slice[n],
                      #  sub_slice[n+1]
                        ) 
                
                to_remove.extend(tuple_)
                tuple_bond_order_ports.append(tuple_)
                
            # Remove stored elements from the list
            for element in to_remove:
                if element in PLP[Chain_]:  # Check before removing
                    PLP[Chain_].remove(element)
                elif element in PLP[Chain_ + 1]:
                    PLP[Chain_ + 1].remove(element)
                    
                
        
    Ports_list_crosslinker_tmp = []
    for n in range(total_crosslinkers_Degree_Crosslinking):
        each_crosslinker_ports = []
        for p in  range(int(crosslink_sites_per_molecule * 2)):   # 2 crosslinking sites per side of crosslinker * 2 
            name_port_it = f"{Crosslinker_port_name}-{n}_{p}"
            each_crosslinker_ports.append(name_port_it)
            
            
        Ports_list_crosslinker_tmp.append(each_crosslinker_ports)
        
    tuple_bond_order_sites = [tuple(sublist) for sublist in Ports_list_crosslinker_tmp]
    
    
    print(tuple_bond_order_sites)
    print(tuple_bond_order_ports)
    if len(tuple_bond_order_sites) != len(tuple_bond_order_ports):
        raise ValueError("The number of sites tuples does not match the number of ports tuples")

        

    return tuple_bond_order_ports, tuple_bond_order_sites




def Crosslinker_Child_Remove_atoms(Polymer_compound = mb.Compound(),
                               child_name = "CRS-0",
                               anchor_atom = "C",
                               anchor_atom_n_bonds = 3,
                               target_atom = ["H", "O"],
                               atoms_to_prune = 2,
                               port_name = "Port",
                               bond_length_factor = 2,
                               prune_extra = False,
                               list_extra_pruning = ["O",["H"]],
                               list_to_exclude = [],
                              
                               ):
    
    c = 0
    
    #######
    print("Non nested crosslinker")
    #######
    removed_particle_dictionary =  dict()
    for child in Polymer_compound.children:
        if child.name  == child_name:
            list_alpha = []
            for particle in child.particles():
                if c < atoms_to_prune:
                    #print(f"Found the {anchor_atom} atom")
                    proceed = False
                    if particle.name ==anchor_atom:
                        list_names = []
                        list_positions = []
                        list_id = []
                        particle_list = []
                      
                        proceed = True
                                            
                        for p in particle.direct_bonds():
                            list_names.append(p.name)
                            list_positions.append(p.pos)
                            list_id.append(p)
                            particle_list.append(p)
                            
                        to_remove = []
                        if (len(list_names) == anchor_atom_n_bonds) and all(elem in list_names for elem in target_atom) and proceed ==True :
                            for Ta in target_atom:
                                if prune_extra == True:
                                    for particle_attached in particle_list:
                                        #print("particles attached")
                                        #print(particle_attached)
                                        if particle_attached.name == list_extra_pruning[0]:
                                            for particle_attached_target in particle_attached.direct_bonds():
                                                for element_to_prune_attached_target in list_extra_pruning[1]:
                                                    if particle_attached_target.name == element_to_prune_attached_target:
                                                       # print(f"Removing {particle_attached_target.name} from the {particle_attached.name }")
                                                        to_remove.append(particle_attached_target)
                                else:
                                    pass
                                    
                                position_H = list_positions[list_names.index(Ta)]
                                # Define the positions of the two points
                                point1 = position_H
                                point2 = particle.pos
                                if Ta not in list_to_exclude:
                                    # Calculate the vector from point1 to point2
                                    vector = point1 - point2
                                    # Calculate the magnitude (length) of the vector
                                    magnitude = np.linalg.norm(vector)
                                    # Calculate the unit vector (normalize the vector)
                                    unit_vector = vector / magnitude
                                    distance = np.linalg.norm(point2 - point1)
                                    name_ = port_name +"_"+ str(c)
                                    print(f"added port : {name_}")
                                    Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                    to_remove.append(list_id[list_names.index(Ta)])
                                    #Polymer_compound.remove(list_id[list_names.index(Ta)])
                                    Polymer_compound.remove(to_remove)
                                    list_alpha.append((c,point2,point1))
                                    c+=1
                                else:
                                    # Calculate the vector from point1 to point2
                                    vector = point1 - point2
                                    # Calculate the magnitude (length) of the vector
                                    magnitude = np.linalg.norm(vector)
                                    # Calculate the unit vector (normalize the vector)
                                    unit_vector = vector / magnitude
                                    distance = np.linalg.norm(point2 - point1)
                                    ######
                                    ###
                                    #
                                    # Rotation angle (45 degrees)
                                    theta = np.radians(45)

                                    # Rotation matrix for XY plane (Z-axis rotation)
                                    R_XY = np.array([
                                        [np.cos(theta), -np.sin(theta), 0],
                                        [np.sin(theta),  np.cos(theta), 0],
                                        [0,             0,             1]
                                    ])

                                    # Rotation matrix for XZ plane (Y-axis rotation)
                                    R_XZ = np.array([
                                        [np.cos(theta), 0, np.sin(theta)],
                                        [0,             1, 0],
                                        [-np.sin(theta), 0, np.cos(theta)]
                                    ])

                                    # Apply the rotations
                                    vector_rotated_xy = R_XY @ unit_vector  # Rotate in XY plane
                                    vector_rotated_xz = R_XZ @ vector_rotated_xy  # Rotate the result in XZ plane
                                    #
                                    ###
                                    #####
                                    name_ = port_name +"_"+ str(c)
                                    Polymer_compound.add(mb.Port(anchor=particle, orientation=vector_rotated_xz, separation=distance * bond_length_factor), name_)
                                    #to_remove.append(list_id[list_names.index(Ta)])
                                    #Polymer_compound.remove(list_id[list_names.index(Ta)])
                                    #Polymer_compound.remove(to_remove)
                                    list_alpha.append((c,point2,point1))
                                    c+=1
                                    
                             
                                #print("Removed particle")
                            
                                                            
                else:
                    break
                
    removed_particle_dictionary[child_name] = list_alpha
            

                
    return Polymer_compound, removed_particle_dictionary


def Crosslinker_Child_Remove_atoms_NESTED(Polymer_compound = mb.Compound(),
                                backbone_T = dict(),
                               child_name = "CRS-0",
                               atoms_to_prune = 2,
                               port_name = "Port",
                               bond_length_factor = 2,
                               
                               ):
    
    c = 0
    
    #######
    #print("Doing crosslinkers")
    removed_particle_dictionary =  dict()
    for child in Polymer_compound.children:
        for nested_dictionary in backbone_T:
            #print(f"doing dictionary : {nested_dictionary}")
            anchor_atom = nested_dictionary["anchor_atom"]
            anchor_atom_n_bonds = nested_dictionary["anchor_atom_n_bonds"]
            target_atom = nested_dictionary["target_atom"]
            list_extra_pruning = nested_dictionary["extra_pruning"]
            if list_extra_pruning:
                prune_extra = True
            else:
                prune_extra = False
            list_to_exclude = nested_dictionary["exclude"]

            if child.name  == child_name:
                list_alpha = []
                for particle in child.particles():
                    if c < atoms_to_prune:
                        #print(f"Found the {anchor_atom} atom")
                        proceed = False
                        if particle.name ==anchor_atom:
                            list_names = []
                            list_positions = []
                            list_id = []
                            particle_list = []
                        
                            proceed = True
                                                
                            for p in particle.direct_bonds():
                                list_names.append(p.name)
                                list_positions.append(p.pos)
                                list_id.append(p)
                                particle_list.append(p)
                                
                            to_remove = []
                            if (len(list_names) == anchor_atom_n_bonds) and all(elem in list_names for elem in target_atom) and proceed ==True :
                                for Ta in target_atom:
                                    if prune_extra == True:
                                        for particle_attached in particle_list:
                                            #print("particles attached")
                                            #print(particle_attached)
                                            if particle_attached.name == list_extra_pruning[0]:
                                                for particle_attached_target in particle_attached.direct_bonds():
                                                    for element_to_prune_attached_target in list_extra_pruning[1]:
                                                        if particle_attached_target.name == element_to_prune_attached_target:
                                                        # print(f"Removing {particle_attached_target.name} from the {particle_attached.name }")
                                                            to_remove.append(particle_attached_target)
                                    else:
                                        pass
                                        
                                    position_H = list_positions[list_names.index(Ta)]
                                    # Define the positions of the two points
                                    point1 = position_H
                                    point2 = particle.pos
                                    if Ta not in list_to_exclude:
                                        # Calculate the vector from point1 to point2
                                        vector = point1 - point2
                                        # Calculate the magnitude (length) of the vector
                                        magnitude = np.linalg.norm(vector)
                                        # Calculate the unit vector (normalize the vector)
                                        unit_vector = vector / magnitude
                                        distance = np.linalg.norm(point2 - point1)
                                        name_ = port_name +"_"+ str(c)
                                        #print(f"added port : {name_}")
                                        Polymer_compound.add(mb.Port(anchor=particle, orientation=unit_vector, separation=distance * bond_length_factor), name_)
                                        to_remove.append(list_id[list_names.index(Ta)])
                                        #Polymer_compound.remove(list_id[list_names.index(Ta)])
                                        Polymer_compound.remove(to_remove)
                                        list_alpha.append((c,point2,point1))
                                        c+=1
                                    else:
                                        # Calculate the vector from point1 to point2
                                        vector = point1 - point2
                                        # Calculate the magnitude (length) of the vector
                                        magnitude = np.linalg.norm(vector)
                                        # Calculate the unit vector (normalize the vector)
                                        unit_vector = vector / magnitude
                                        distance = np.linalg.norm(point2 - point1)
                                        ######
                                        ###
                                        #
                                        # Rotation angle (45 degrees)
                                        theta = np.radians(45)

                                        # Rotation matrix for XY plane (Z-axis rotation)
                                        R_XY = np.array([
                                            [np.cos(theta), -np.sin(theta), 0],
                                            [np.sin(theta),  np.cos(theta), 0],
                                            [0,             0,             1]
                                        ])

                                        # Rotation matrix for XZ plane (Y-axis rotation)
                                        R_XZ = np.array([
                                            [np.cos(theta), 0, np.sin(theta)],
                                            [0,             1, 0],
                                            [-np.sin(theta), 0, np.cos(theta)]
                                        ])

                                        # Apply the rotations
                                        vector_rotated_xy = R_XY @ unit_vector  # Rotate in XY plane
                                        vector_rotated_xz = R_XZ @ vector_rotated_xy  # Rotate the result in XZ plane
                                        #
                                        ###
                                        #####
                                        name_ = port_name +"_"+ str(c)
                                        Polymer_compound.add(mb.Port(anchor=particle, orientation=vector_rotated_xz, separation=distance * bond_length_factor), name_)
                                        #to_remove.append(list_id[list_names.index(Ta)])
                                        #Polymer_compound.remove(list_id[list_names.index(Ta)])
                                        #Polymer_compound.remove(to_remove)
                                        list_alpha.append((c,point2,point1))
                                        c+=1
                                        
                                
                                    #print("Removed particle")
                                
                                                                
                    else:
                        break
                
    removed_particle_dictionary[child_name] = list_alpha
            

                
    return Polymer_compound, removed_particle_dictionary

########################################## Crosslinking degree #########################################

def round_to_next_multiple_of_4(n):
    if n > 4:
        if n % 4 == 0:
            return n, True
        else:
            return ((n // 4) + 1) * 4, False
        
    return None, None  # Return None if n is 4 or less

def determine_max_crosslinker(ru = int, sites = int,  dividend = int, NumChains = int):
    val = ((ru*sites)/4) * (NumChains-1)
    return val

def is_odd(num):
    return isinstance(num, int) and num % 2 == 1
 

def determine_number_crosslinks( repeat_units=6,  number_of_polymers=3,  sites_per_repeat_unit=1,  crosslinks_unit_per_chain=4, straight_to_crosslink_percentage = True, Num_Crosslink_Compounds_between_Chains = 1, mute = True):
    
    if number_of_polymers <= 1 or repeat_units < 4:
        raise ValueError("Number of polymers must be greater than 1, number of repeat units must be greater than 4")

    else:
        total_sites = repeat_units * sites_per_repeat_unit * number_of_polymers  
        sites_available_for_crosslinker = int(crosslinks_unit_per_chain/2)
        sites_per_chain = (repeat_units * sites_per_repeat_unit)//(crosslinks_unit_per_chain/2)
        
        if not mute:
            print(f"You have {number_of_polymers} polymer chains \n each chain has {repeat_units} repeat units \n each repeat unit has {sites_per_repeat_unit} functional sites \n with a total of {total_sites} sites for all the chains")
            
            print(f"The crosslinker takes {sites_available_for_crosslinker} sites per chain")
            
            print(f"Each chain then has {sites_per_chain} crosslinking sites")
        
        if straight_to_crosslink_percentage == False:
                    
            if number_of_polymers == 2:
                
                max_number_of_crosslinks = (sites_per_repeat_unit * repeat_units) 
                max_per_top_bottom = max_number_of_crosslinks//(crosslinks_unit_per_chain/2)
                max_per_mid = int(max_number_of_crosslinks//(crosslinks_unit_per_chain/2))
                
            elif number_of_polymers > 2:
                
                if crosslinks_unit_per_chain ==4:
                    Repeat_unit_calc, exact_boolean = round_to_next_multiple_of_4(repeat_units)
                else:
                     Repeat_unit_calc = repeat_units
                     exact_boolean = True
                
                if exact_boolean == True:
                    
                    max_number_of_crosslinks = int(((sites_per_repeat_unit * Repeat_unit_calc) /crosslinks_unit_per_chain) * (number_of_polymers-1)) 
                    
                    max_per_top_bottom = (Repeat_unit_calc//crosslinks_unit_per_chain)
                    max_per_mid = ((Repeat_unit_calc//crosslinks_unit_per_chain)) 
                    
                elif exact_boolean == False:
                    
                    if not mute:
                        print("Not an exact division between RU/Cross")
                    
                    max_number_of_crosslinks = int(((sites_per_repeat_unit * Repeat_unit_calc) /crosslinks_unit_per_chain) * (number_of_polymers-1)) -1
                    max_per_top_bottom = (Repeat_unit_calc//crosslinks_unit_per_chain)
                    max_per_mid = ((Repeat_unit_calc//crosslinks_unit_per_chain)) - 1
                    
        else:
            
            max_number_of_crosslinks = Num_Crosslink_Compounds_between_Chains * (number_of_polymers-1)
            max_per_top_bottom = Num_Crosslink_Compounds_between_Chains//2
            max_per_mid = Num_Crosslink_Compounds_between_Chains
            
                
        used_fg = int(max_number_of_crosslinks * crosslinks_unit_per_chain)
        percentage_CR = int((used_fg/ total_sites) *100)
        
        if straight_to_crosslink_percentage == False:
            if not mute:
                print("For a system with all polymer chains connected with the maximum number of crosslink agents in an x,y grid")
                print(f"    The maximum number of crosslink molecules to be used is: {max_number_of_crosslinks} ")
                print(f"    Which accounts for {used_fg} functional sites used")
                print(f"    That is {percentage_CR} % crosslinking")
                
                print(f"With bonds between polymer chains as:")
        
        else:
            if not mute:
                print("For a system with all polymer chains connected with the maximum number of crosslink agents in an x,y grid")
                print(f"    The {Num_Crosslink_Compounds_between_Chains} crosslinking compounds between chains \n accounts for  {percentage_CR} % crosslinking")
                
                print(f"With bonds between polymer chains as:")
            
            
        num_bonds_between_Chains_tuple = (max_per_top_bottom, max_per_mid)
        dict_bonds_CRS = {}
        my_list = list(range(number_of_polymers-1)) 
        for n in my_list:
            if n == 0:  
                dict_bonds_CRS[f"POL-{n}/POL-{n+1}"] = max_per_top_bottom
                
            elif n != 0 and n != my_list[-1]:
           
                dict_bonds_CRS[f"POL-{n}/POL-{n+1}"] = max_per_mid
                
            elif   n == my_list[-1]:
                
                if is_odd(number_of_polymers) == True:
            
                    dict_bonds_CRS[f"POL-{n}/POL-{n+1}"] = max_per_mid
                
                else:
                    dict_bonds_CRS[f"POL-{n}/POL-{n+1}"] = max_per_top_bottom
                    

        if not mute:   
            print(dict_bonds_CRS)

                
    return max_number_of_crosslinks, percentage_CR, dict_bonds_CRS, num_bonds_between_Chains_tuple 

 
def get_closest_crosslinker(crosslink_dict, target_percentage):
    closest_key = min(crosslink_dict, key=lambda k: abs(crosslink_dict[k] - target_percentage))
    closest_value = crosslink_dict[closest_key]
    return closest_key, closest_value


def CRP_num_between_chains(repeat_units=6,  number_of_polymers=3,  sites_per_repeat_unit=1,  crosslinks_unit_per_chain=4, desired_crosslinking_percentage = 1, mute = True ):
    
    max_number_of_crosslinks, _, _, num_bonds_between_Chains_tuple = determine_number_crosslinks( repeat_units=repeat_units,
                                                                    number_of_polymers=number_of_polymers,
                                                                    sites_per_repeat_unit=sites_per_repeat_unit,
                                                                    crosslinks_unit_per_chain=crosslinks_unit_per_chain,
                                                                    straight_to_crosslink_percentage = False,
                                                                    Num_Crosslink_Compounds_between_Chains =1, 
                                                                    mute = mute,
                                                                    )
    dict_crosslinking = {}
        
    for i in list(range(1, num_bonds_between_Chains_tuple[1] + 1)):

        _, percentage_CR_1, _, _ = determine_number_crosslinks( repeat_units=repeat_units,
                                                                number_of_polymers=number_of_polymers,
                                                                sites_per_repeat_unit=sites_per_repeat_unit,
                                                                crosslinks_unit_per_chain=crosslinks_unit_per_chain,
                                                                straight_to_crosslink_percentage = True,
                                                                Num_Crosslink_Compounds_between_Chains = i,
                                                                mute = mute,
                                                                )
        
        dict_crosslinking[i] =   percentage_CR_1
        
        
    closest_monomer_match, actual_match = get_closest_crosslinker(dict_crosslinking,desired_crosslinking_percentage )
    
    if number_of_polymers <=2:
        total_monomer = closest_monomer_match 
    else:
        total_monomer = closest_monomer_match * (number_of_polymers-1)
        
        

    keys_monomers = list(dict_crosslinking.keys())
    if keys_monomers and keys_monomers[-1] == closest_monomer_match:
        boolean_highest_degree_CR_possible = True
    else:
        boolean_highest_degree_CR_possible = False
        
        
    return dict_crosslinking , closest_monomer_match , total_monomer, actual_match, boolean_highest_degree_CR_possible

######################
###########
###
def Crosslink_Pipeline(
    main_name = "ethanol",                                   # Monomer to polimerize gives PVA ; not vinylic alcohol as it's insaturated
    MB_molecule_to_create_monomer = Alcohol(chain_length=2), # Mbuild compound, ethanol molecule
    fetch_compound_from_database = False,
    first_port_monomer = "COH",
    second_port_monomer = "CH3",
    
    repeat_units = 100,                                       # Repeat Units in polymer chain
    number_of_polymers = 3,                                   # Number of polymer chains inside the box
    depth_value = 1,                                          # Spacing between polymer chains in box
    Polymer_residue_name = "POL",                             # name of polymer chain residue
    Polymer_port_name = "Port",                               # name of polymer chain backbone ports
    Dummy_atoms = False,
    Dummy_atom_dictionary = None,

                                                              #Crosslinker details
    crosslinker_name = "Glutaraldehyde",                      # Name of crosslinker agent
    Crosslinker_residue_name = "CRS",                         # name of crosslinker residue name
    Crosslinker_port_name = "Site",                           # name of cr chain backbone ports
    Crosslinker_compound = None,                              # Compound (mb.Compound)


                                                              #Crosslinking details
    sites_per_repeat_unit = 1,                                # Number of functional groups (f.ex hydroxyls in PVA in each chain)
    desired_CR_degree = 50,                                   # Degree in %, that percentage (f.ex number of hydroxyls attached to the crosslinker in PVA, not hydroxyls anymore)
    
    
    polymer_backbone_target_list = list(),
    
    crosslinker_backbone_target_list = list(),
    
    prune_extra_flag =False,                                     # No pruning extra atoms
    crosslink_sites_per_molecule = 1,
    
    Heterogeneous_crosslinking = False,
    
    dimension_0_1_2 = 0,
    
    
):
    Randomness = True  # Boolean variable, meaning the ports in the polymer chain are going to be selected randomly, not continously
    #####################Parameters#################################################
    #MB_molecule_to_create_monomer = Alcohol(chain_length=2)
    #MB_molecule_to_create_monomer.visualize()

    tuple_sites = [
        
        ("C=O", "site_1", 2), 
        ("C=O", "site_3", 2), 
        ("CO_g", "site_2", 2), 
        ("CO_g", "site_4", 2), 
        
    ]

   
    


    ############################################
    ################### Code Details ###########
    '''
    Polymer Chain

    Take the molecule and select target ports for the polymer_builder function in mbuild to use this molecule, create a monomer and create a polymer chain
        1. The molecule can be fetched from the PubChem database, and loaded as an mbuild compound:
            fetch_compound_from_database = True
            main_name = "Name of molecule"
            MB_molecule_to_create_monomer = None
            
        2. The molecule can be passed to the function in the form of an mbuild monomer: 
            fetch_compound_from_database = False
            main_name = "Name of molecule"
            MB_molecule_to_create_monomer = mb.Compound()
        
    Create a monomer that can be polymerized
    Example:
        PVA from ethanol;
        1. An H atom has to be removed from the C(C)(H)(H)(O)H -> Port "COH" to create a port and remove that hydrogen atom; first_port_monomer
        2. An H atom has to be removed from the C(C)(H)(H)(H) -> Port "CH3" to create a port and remove that hydrogen atom; second_port_monomer


    General 

    1. Each polymer chain is added to an MB Box with the residue name of the chain being:
        POL-0, POL-1, POL-3 ... where the number corresponds to the polymer chain index, and the Polymer_residue_name
        
    2. The Hydrogen atoms in the hydroxyl group of the PVA molecule are pruned, the resulting port is given the name:
        Port-0-0, Port-0-1... where the first number corresponds to the polymer chain index, and the latter to the port number (oxygen atom in the polymer chain that was a hydroxyl) : Polymer_port_name
        
    3. The crosslinker is created, and the ports correspond to: 
        Site_1 and Site_2 in one side, and Site_3 and Site_4 on the other side. Refer to Parameters tuple_sites
        
    4. The  crosslink_sites_per_molecule signals the code how many sites per polymer chain backbone the crosslinker attaches to,
        f.ex glutaraldehyde attaches to 2 sites in each polymer chain  (PVA)
        
    Crosslinking details


    '''
    #####################################################################################################################################################################################################################
    dict_crosslinkers_to_percentage, crosslinkers_per_pair_of_chains, total_crosslinkers_Degree_Crosslinking, ADC, boolean_highest_degree_CR_possible  = CRP_num_between_chains( repeat_units=repeat_units,
                                                                                                                                        number_of_polymers=number_of_polymers,
                                                                                                                                        sites_per_repeat_unit= sites_per_repeat_unit,
                                                                                                                                        crosslinks_unit_per_chain= crosslink_sites_per_molecule*2, # It's two per chain, 4 for a complete bridge between two chains
                                                                                                                                        desired_crosslinking_percentage = desired_CR_degree, 
                                                                                                                                    
                                                                                                                                        )
    '''
    This lines determines the number of ports to add to each polymer chain based on how many ports the crossliner needs to attach;
        f.ex PVA-Glutaraldehyde
            Each crosslinker molecule needs two Hydroxyl oxygens in each Polymer Chain to attach,
            so if you need 3 crosslinker molecules between each pair of chains to attain the desire degree of crosslinking
            for each crosslinker molecule you need 2 ports in each chain
    '''
    if Heterogeneous_crosslinking:
        tuple_case_ = False
       
         
        max_crosslinks = ((number_of_polymers - 1) * repeat_units * sites_per_repeat_unit)
        desired_crosslinks = round(max_crosslinks * (desired_CR_degree / 100))
        max_crosslinks_between_chains = int(max_crosslinks/(number_of_polymers - 1))
        
        total_sites = repeat_units * sites_per_repeat_unit * number_of_polymers
        sites_p_chain = repeat_units * sites_per_repeat_unit
        
        crosslinkers_per_pair = round(desired_crosslinks / (number_of_polymers - 1))
        Number_of_ports_per_chain = crosslinkers_per_pair
     
        
        print(f"You have a total of {total_sites} sites")
        print(f"Each chain has {sites_p_chain} sites")
        print(f"For 100% degree of crosslinking you need {max_crosslinks} crosslinkers")
        print(f"Which means {max_crosslinks_between_chains} between chains")
        #print(f"For a degree of {desired_CR_degree}% you need {desired_crosslinks} crosslinkers")
        #print(f"Number of ports per chain: {Number_of_ports_per_chain}\n")


        dict_Het = {}
        for i in list(range(1, max_crosslinks_between_chains + 1)):
            
            dict_Het[i] = 100* (i* (number_of_polymers)/total_sites)
            
        for key, val in dict_Het.items():
            print(f"For {key} crosslinker/s between chains you get {val}% crosslinking")
        
            
            
        
        if boolean_highest_degree_CR_possible:
            print(f'Which is the highest crosslinking percentage attainable')
            
        
        # Find the key with the closest value
        closest_key = min(dict_Het, key=lambda k: abs(dict_Het[k] - desired_CR_degree))
        closest_value = dict_Het[closest_key]
        total_crosslinkers_Degree_Crosslinking = int((number_of_polymers - 1) * closest_key)
         
        print(f"The closest value for the desired CR degreee ({desired_CR_degree}%) \n is {closest_value}, \nwhich accounts for {closest_key} crosslinker between chains \n with a total of {total_crosslinkers_Degree_Crosslinking} crosslinkers")
        
        Number_of_ports_per_chain = closest_key
        
        total_crosslinkers_Degree_Crosslinking = (int((number_of_polymers - 1) * closest_key)) 
        
        crosslinkers_per_pair_of_chains = closest_key
        
    else:
        Number_of_ports_per_chain = crosslinkers_per_pair_of_chains * crosslink_sites_per_molecule  
        if crosslink_sites_per_molecule ==2:
            tuple_case_ = True
        elif crosslink_sites_per_molecule ==1:
            tuple_case_ = False

        print(f"Number of ports per chain {Number_of_ports_per_chain}")
        ################################################# Print Statement #########################################################################################################################################################################################################
        print("This values are the number of crosslinking monomers that will be bonded between each pair of chains to attain the desire degree of crosslinking")
        print(dict_crosslinkers_to_percentage)
        print(f"To attain the degree of crosslinking {desired_CR_degree}% \n {total_crosslinkers_Degree_Crosslinking} crosslinkers are needed \n {crosslinkers_per_pair_of_chains} between each polymer chain")
        print(f"The actual attainable value of crosslinking is {ADC}%")
        if boolean_highest_degree_CR_possible:
            print(f'Which is the highest crosslinking percentage attainable')
    #####################################################################################
    #############################################
    ###################################
    main_polymer_bead, end_bead = Builder_pipeline_Polymer(main_name,
                                                           first_port = first_port_monomer,
                                                           second_port = second_port_monomer,
                                                           add_ports = True,
                                                           rigid_port = False,
                                                           fetch_compound = fetch_compound_from_database,
                                                           MB_compound = MB_molecule_to_create_monomer,
                                                           factor_length = 1.2)
    
    polymer_compound_ = mb.recipes.Polymer(monomers=[main_polymer_bead])
    polymer_compound_.build( n=repeat_units,)
    polymer_compound = mb.Compound()
    polymer_compound.add(polymer_compound_)
    polymer_compound.visualize()
    
    
                    
    polymer_compound.visualize(show_ports = True)
    print("Done with loading compounds")
    if Crosslinker_compound:
        base_crosslinker = Crosslinker_compound
    else:
        digits, cid = fetch_and_plug_CR(crosslinker_name)
        smiles, mol = get_mol_CR(cid)
        base_crosslinker = Build_3D_CR(smiles) 

    #base_crosslinker.visualize(show_ports = True) 
    Compound_0  = mb.Compound()
    list_of_polymers = []
    polymer_target_atom_positions = {}
    for i in range(number_of_polymers):
            polymer = mb.clone(polymer_compound)
            polymer.name = f"{Polymer_residue_name}-{i}"
            polymer_target_atom_positions[f"{Polymer_residue_name}-{i}"] = []
            list_of_polymers.append(polymer)
            Compound_0.add(polymer, label = f"{Polymer_residue_name}-{i}")

    list_of_crosslinkers = []
    for i in range(total_crosslinkers_Degree_Crosslinking):
        crs = mb.clone(base_crosslinker)
        crs.name = f"{Crosslinker_residue_name}-{i}"
        list_of_crosslinkers.append(crs)
        Compound_0.add(crs, label = f"{Crosslinker_residue_name}-{i}")



    if Heterogeneous_crosslinking:
        shift_value = Find_Shift_value (
                            Monomer = main_polymer_bead,
                            targets_dictionary_list = polymer_backbone_target_list,
                            dimension = dimension_0_1_2,
        
                            )
    else:
        shift_value = 0
        
    for i,child in enumerate(Compound_0.children):
        if f"{Crosslinker_residue_name}" not in child.name:
            if child.name == f"{Polymer_residue_name}-0":
                first_positions = child.xyz
            

            if child.name != f"{Polymer_residue_name}-0":
                new_array = first_positions
                new_array[:,2] = first_positions[:,2]+ (depth_value)
                if Heterogeneous_crosslinking == True: # New added flag
                    new_array[:,0] = first_positions[:,0]+ (shift_value)
                    
                child.xyz = new_array
                first_positions = child.xyz
                
               
    polymer_target_atom_positions_final, _  = Polymer_anchor_atoms(Polymer_compound = Compound_0,
                                                                    Target_atom_positions = polymer_target_atom_positions,
                                                                    polymer_backbone_targets = polymer_backbone_target_list,
                                                                    bonded_neighbors_bool = True,
                                                        
                                                                    )
    print(polymer_target_atom_positions_final) 
    print(polymer_target_atom_positions_final["POL-0"])                                                           
    print("Done getting atom positions")
    print(polymer_target_atom_positions_final)
    print("Those are the positions")
    if len(list_of_polymers)<=2:
        print("Two polymer chains")
        for i in range(len(list_of_polymers)):
        
            if i ==0:
                groups_pol_1 = select_continuous_groups_random(polymer_target_atom_positions_final,
                                                               f"{Polymer_residue_name}-0",
                                                               Number_of_ports_per_chain,
                                                               crosslink_sites_per_molecule,
                                                               number_of_polymers,
                                                               force_even_start=boolean_highest_degree_CR_possible,
                                                               tuple_case=tuple_case_)
                keys_list = list(groups_pol_1.keys())
                sorted_dict = dict(sorted(groups_pol_1.items()))
                groups_pol = sorted_dict
                factor =1
                key_list = keys_list
                
            elif  i == len(list_of_polymers) - 1:
                factor = 1
                groups_pol = extract_matching_tuples(polymer_target_atom_positions_final, f"{Polymer_residue_name}-{i}") ## changed key_list to  New_key_list
            
            num_keys = len(groups_pol)
            compound_iterations, removed_particle_dictionary  = Polymer_Child_Remove_atoms(Polymer_compound =Compound_0,
                                                                                    child_name = f"{Polymer_residue_name}-{i}",
                                                                                    anchor_atom = polymer_backbone_target['anchor_atom'],
                                                                                    anchor_atom_n_bonds = polymer_backbone_target['anchor_atom_n_bonds'],
                                                                                    target_atom = polymer_backbone_target['target_atom'],
                                                                                    atoms_to_prune = Number_of_ports_per_chain*factor,
                                                                                    port_name = f"{Polymer_port_name}-{i}",
                                                                                    bond_length_factor = 1, 
                                                                                    dictionary_port_index = groups_pol,
                                                                        )
        

    else:
        lap = []
        print("More than two polymer chains")
        for i in range(len(list_of_polymers)):
            print(f"Chain number {i}")
            if i ==0:
                if Heterogeneous_crosslinking:
                    factor =1
                else: 
                    factor = 1
                Flag = False
                attached_ports = {}
                #######################################################
                
                if Heterogeneous_crosslinking:
                    # Extract the list of tuples for the selected key
                    selected_data = polymer_target_atom_positions_final[f"{Polymer_residue_name}-0"]

                    # Convert list of tuples into a new dictionary
                    groups_pol_1 = {idx: pos for idx, pos in selected_data}
                    sorted_dict = dict(sorted(groups_pol_1.items(), key=lambda item: item[1][0], reverse=True))
                    groups_pol = sorted_dict
                    
                else:
                    #Number_of_ports_per_chain = This is the number of ports to vacate per chain, so if 2 monomers to add between chains then that's 2 * crosslink_sites_per_molecule (2) which gives 4
                    # crosslink_sites_per_molecule = 2, for the case of Glutaraldehyde how many OH's per each crossliker monomer for individual chain
                    groups_pol_1 = select_continuous_groups_random(polymer_target_atom_positions_final,
                                                                f"{Polymer_residue_name}-0",
                                                                Number_of_ports_per_chain,
                                                                crosslink_sites_per_molecule,
                                                                number_of_polymers,
                                                                    force_even_start=boolean_highest_degree_CR_possible,
                                                                tuple_case=tuple_case_,
                                                                )
                
                    sorted_dict = dict(sorted(groups_pol_1.items(), key=lambda item: item[1][0], reverse=True))
                    groups_pol = sorted_dict
                
                print("First groups selected")
                print(len(groups_pol))
                print(groups_pol)
                #print(f"{i+1} Chain")
                #print(f"Keys list : {keys_list}")
                #print(groups_pol)
                
            elif i> 0 and i< len(list_of_polymers) - 1: 
                factor = 2
                Flag = True
                if Heterogeneous_crosslinking:

                    filtered_dict_F = extract_matching_tuples(polymer_target_atom_positions_final, f"{Polymer_residue_name}-{i}") ## changed key_list to  New_key_list 
                    #removed_keys, selected_pairs, selected_pairs_dict = process_dicts(previous_dict, groups_pol_1, tuple_case = tuple_case_)
                    #filtered_dict = filter_selected_pairs(selected_pairs_dict, attached_ports) # Remove the ports that are in front of the new chain, as they are already going to be added
                    #filtered_dict_F  = filter_non_overlapping_continuous_pairs(filtered_dict, tuple_case = tuple_case_)
                    
                    
                else:
                    groups_pol_1 = extract_matching_tuples(polymer_target_atom_positions_final, f"{Polymer_residue_name}-{i}") ## changed key_list to  New_key_list 
                    removed_keys, selected_pairs, selected_pairs_dict = process_dicts(previous_dict, groups_pol_1, tuple_case = tuple_case_)
                    filtered_dict = filter_selected_pairs(selected_pairs_dict, attached_ports) # Remove the ports that are in front of the new chain, as they are already going to be added
                    filtered_dict_F  = filter_non_overlapping_continuous_pairs(filtered_dict, tuple_case = tuple_case_)
                    
                    #print(f"{i+1} Chain")
                    #print(f"Keys list : {New_keys_list}")
                    #print(f"New numbers {filtered_dict}")
                
                
                groups_pol = filtered_dict_F
                
                #print("Filtered Dict")
                #print(filtered_dict)
                #print(groups_pol)
                
                
            elif  i == len(list_of_polymers) - 1:
                factor = 1
                Flag = True
                
                if Heterogeneous_crosslinking:
                    groups_pol = extract_matching_tuples(polymer_target_atom_positions_final, f"{Polymer_residue_name}-{i}") ## changed key_list to  New_key_list 
                     
                else:
                    groups_pol_1 = extract_matching_tuples(polymer_target_atom_positions_final, f"{Polymer_residue_name}-{i}") ## changed key_list to  New_key_list 
                    removed_keys, selected_pairs, selected_pairs_dict = process_dicts(previous_dict, groups_pol_1, tuple_case = tuple_case_)
                    #filtered_dict_F = {k: v for k, v in selected_pairs_dict.items() if k + 1 in selected_pairs_dict or k - 1 in selected_pairs_dict} 
                    #filtered_dict = filter_selected_pairs(selected_pairs_dict, attached_ports) # Remove the ports that are in front of the new chain, as they are already going to be added
                    #print(f"{i+1} Chain, Last")
                    #print(f"Keys list : {new_numbers}") 
                    groups_pol = selected_pairs_dict
                
            
                
            
            
            if isinstance(polymer_backbone_target_list, list):
                if len(polymer_backbone_target_list) == 1:
                    # Unwrap the list to get the single dictionary
                    polymer_backbone_target = polymer_backbone_target_list[0]
            
                    compound_iterations, removed_particle_dictionary, attached_ports, L_A_P  = Polymer_Child_Remove_atoms(Polymer_compound =Compound_0,
                                                                                            child_name = f"{Polymer_residue_name}-{i}",
                                                                                            anchor_atom = polymer_backbone_target['anchor_atom'],
                                                                                            anchor_atom_n_bonds = polymer_backbone_target['anchor_atom_n_bonds'],
                                                                                            target_atom = polymer_backbone_target['target_atom'],
                                                                                            atoms_to_prune = Number_of_ports_per_chain*factor,
                                                                                            port_name = f"{Polymer_port_name}-{i}",
                                                                                            bond_length_factor = 1, 
                                                                                            dictionary_port_index = groups_pol,
                                                                                            restricted_insertion = Flag, 
                                                                                            ports_dict = attached_ports,
                                                                                            
                                                                                )
                    
                elif len(polymer_backbone_target_list) == 2 :
                    
                    print(f"Iteration {i}")
                    
                    compound_iterations, removed_particle_dictionary, attached_ports, L_A_P  = Polymer_Child_Remove_atoms_NESTED(Polymer_compound =Compound_0,
                                                                                            child_name = f"{Polymer_residue_name}-{i}",
                                                                                            targets_dictionary_list= polymer_backbone_target_list,
                                                                                                
                                                                                            atoms_to_prune = Number_of_ports_per_chain *factor ,
                                                                                            port_name = f"{Polymer_port_name}-{i}",
                                                                                            bond_length_factor = 1, 
                                                                                            dictionary_port_index = groups_pol,
                                                                                            restricted_insertion = Flag, 
                                                                                            ports_dict = attached_ports,
                                                                                          
                                                                                            
                                                                                )
                    print("attached ports"
                    )
                    print(attached_ports)
                else:
                    print("error polymer")
                    raise TypeError("crosslinker_backbone_target must be a dict or list of dicts")
                
            lap.append(L_A_P)
            previous_dict = groups_pol
            
            
            
        
                                
    print("Done with polymer")
    for i in range(len(list_of_crosslinkers)):

        if isinstance(crosslinker_backbone_target_list, list):
            if len(crosslinker_backbone_target_list) == 1:
                # Unwrap the list to get the single dictionary
                crosslinker_backbone_target = crosslinker_backbone_target_list[0]
                
                compound_iterations, removed_particle_dictionary = Crosslinker_Child_Remove_atoms(Polymer_compound = Compound_0,
                                                                                            child_name = f"{Crosslinker_residue_name}-{i}",
                                                                                            anchor_atom = crosslinker_backbone_target['anchor_atom'],
                                                                                            anchor_atom_n_bonds = crosslinker_backbone_target['anchor_atom_n_bonds'],
                                                                                            target_atom = crosslinker_backbone_target['target_atom'],
                                                                                            atoms_to_prune = crosslink_sites_per_molecule*2,
                                                                                            port_name = f"{Crosslinker_port_name}-{i}",
                                                                                            bond_length_factor = 1,
                                                                                            prune_extra = prune_extra_flag ,
                                                                                            list_extra_pruning = crosslinker_backbone_target["extra_pruning"] , 
                                                                                            list_to_exclude = crosslinker_backbone_target["exclude"]
                                                                                        
                                                                                            )
            elif len(crosslinker_backbone_target_list) == 2 :
                print("multiple sites")
                compound_iterations, removed_particle_dictionary = Crosslinker_Child_Remove_atoms_NESTED(Polymer_compound = Compound_0,
                                                                                                        backbone_T = crosslinker_backbone_target_list,
                                                                                                        child_name = f"{Crosslinker_residue_name}-{i}",
                                                                                                        port_name = f"{Crosslinker_port_name}-{i}",
                                                                                                        bond_length_factor = 1,
                                                                                        
                                                                                        
                                                                                            )
        
        else:
            print("error crosslinker")
            print(crosslinker_backbone_target_list)
            raise TypeError("crosslinker_backbone_target must be a dict or list of dicts")
            
            
    if  tuple_case_ ==True:

        _, tuple_bond_order_sites =  Bond_Order_Tuple(Polymer_port_name, # Port
                                                                        number_of_polymers, # Number of polymer chains 
                                                                        repeat_units, # Number of repeat units in each polymer chain
                                                                        crosslinkers_per_pair_of_chains, # NUmber of crosslinker molecules between each pair of chains
                                                                        crosslink_sites_per_molecule, # Number of functional groups the crosslinker takes from a SINGLE Chain (PVA-Glutarladehyde: 2)
                                                                        total_crosslinkers_Degree_Crosslinking, # Total number of crosslinkers
                                                                        Crosslinker_port_name, # Site
                                                                        )
        #print(tuple_bond_order_sites)
        flattened_lap = [sublist for element in lap for sublist in element]
        #print(flattened_lap)
        #print(flattened_lap)
        grouped_lap = [flattened_lap[i:i+2] for i in range(0, len(flattened_lap), 2)]
        #print(grouped_lap)
        # Process each nested list to extract tuples
        tuple_bond_order_ports = [
            (group[0][i], group[0][i+1], group[1][i], group[1][i+1]) 
            for group in grouped_lap 
            for i in range(0, len(group[0]), 2)
        ]
        #print(tuple_bond_order_ports)
        
    else:
        
        _, tuple_bond_order_sites =  Bond_Order_Non_Tuple(Polymer_port_name, # Port
                                                                        number_of_polymers, # Number of polymer chains 
                                                                        repeat_units, # Number of repeat units in each polymer chain
                                                                        crosslinkers_per_pair_of_chains, # NUmber of crosslinker molecules between each pair of chains
                                                                        crosslink_sites_per_molecule, # Number of functional groups the crosslinker takes from a SINGLE Chain (PVA-Glutarladehyde: 2)
                                                                        total_crosslinkers_Degree_Crosslinking, # Total number of crosslinkers
                                                                        Crosslinker_port_name, # Site
                                                                        )
        #print(tuple_bond_order_sites)
        # Flatten the list
        flattened_lap = [sublist for element in lap for sublist in element]
        #print(flattened_lap)
        grouped_lap = [flattened_lap[i:i+2] for i in range(0, len(flattened_lap), 2)]
        #print(grouped_lap)
        # Process each nested list to extract tuple 
        tuple_bond_order_ports = []
        for sub_list in grouped_lap:
            for ix in range(len(sub_list)-1):
                for iy in range(len(sub_list[ix])):
                    if len(sub_list[ix])>0:
                        tuple_bond_order_ports.append((sub_list[ix][iy], sub_list[ix+1][iy]))
               
            
            
            
            
                       
        
        


    #print(tuple_bond_order_ports)
    ################################### Print Statement ###################################################################
    print(Compound_0)
    #print(tuple_bond_order_ports)
    #print(len(tuple_bond_order_ports))
    #print(tuple_bond_order_sites)
    #print(len(tuple_bond_order_sites))
    ############################################################# Create Crosslink bonds #################################################################
    for i in range(len(tuple_bond_order_ports)):
        #print(f"Crosslinker {Crosslinker_residue_name}-{i}")
        Ports_tuple =  tuple_bond_order_ports[i]
        Sites_tuple =  tuple_bond_order_sites[i]
        
        if len(Ports_tuple) != len(Sites_tuple):
            raise ValueError(f"Mismatch at index {i}: Sites tuple length ({len(Sites_tuple)}) does not match Ports tuple length ({len(Ports_tuple)})")
        
        for c in range(len(Ports_tuple)):        
            #print(f'    {Crosslinker_residue_name}-{i} \n {Sites_tuple[c]} \n {Ports_tuple[c]}')
            
            mb.force_overlap(move_this= Compound_0[f'{Crosslinker_residue_name}-{i}'],
                    from_positions= Compound_0[f'{Sites_tuple[c]}'],
                    to_positions= Compound_0[f'{Ports_tuple[c]}'])
            
    #################################################### Translate compound to positive x, y, z values #####################################################
    translation_vector = 1 - np.min(Compound_0.xyz, axis=0)
    new_pos =  Compound_0.xyz + translation_vector
    Compound_0.xyz  = new_pos
    #################################################### Move Crosslinkers to adequate position #############################################################
    _, list_target_ids  = Polymer_anchor_atoms(Polymer_compound = Compound_0,
                                                Target_atom_positions = polymer_target_atom_positions,
                                                polymer_backbone_targets = polymer_backbone_target_list,
                                                bonded_neighbors_bool= False,
                                                )


    dict_move_CRS = {}

    for child in Compound_0.children:
        for i in range(len(list_of_crosslinkers)):
            if f"{Crosslinker_residue_name}-{i}" == child.name:
                dict_move_CRS[f"{Crosslinker_residue_name}-{i}"] =[]
                for particle in child.particles():
                    for p in particle.direct_bonds():
                        if id(p) in list_target_ids:
                            #print("Value found in list")
                            dict_move_CRS[f"{Crosslinker_residue_name}-{i}"].append(p.pos)
                        else:
                            pass

    averaged_data = {key: np.mean(np.array(value), axis=0) for key, value in dict_move_CRS.items()}
    # Convert dictionary values to a NumPy array
    values = np.array(list(averaged_data.values()))

    # Compute range (max - min) for each dimension
    ranges = np.ptp(values, axis=0)  # ptp gives max - min along axis

    # Identify the dimension with the smallest range
    smallest_change_dimension = np.argmin(ranges)
    
    
    print(averaged_data)
    for child in Compound_0.children:
        for i in range(len(list_of_crosslinkers)):
            #Here I want to add a disclaimer that if i = 0 then I add 0.5 to the dimension with the smallest change (smallest_change_dimension), if it's i=1 then substract 0.5 to the dimension with the smallest change
            # And then add next iteration and then substract. The change should be directly on new_pos
            if f"{Crosslinker_residue_name}-{i}" == child.name:
                #print(f"{Crosslinker_residue_name}-{i} moved to {new_pos}")
                # Apply alternating adjustments
                if i % 2 == 0:
                    adjustment = 0.2
                else:
                     adjustment = -0.2
                    
                new_pos = averaged_data[f"{Crosslinker_residue_name}-{i}"]
                new_pos[smallest_change_dimension] += adjustment
                child.translate_to(new_pos)
                
    Number_Crosslinkers_Used = total_crosslinkers_Degree_Crosslinking
    
    Real_Degree_Crosslinking = ADC
    

    ###########
    if Dummy_atoms == True:
        for atom in Compound_0.particles():
            for key, value in Dummy_atom_dictionary.items():
                if atom.name == key:
                    atom.name = value
                   
                    #print(f"Found dummy atom {key} : replaced it with {value} \n of mass {atom.mass}")
                   
                    
        print("Finished removing dummy atoms")
    ###########
                
    return Compound_0, Number_Crosslinkers_Used, Real_Degree_Crosslinking, Polymer_residue_name, Crosslinker_residue_name,
                
    

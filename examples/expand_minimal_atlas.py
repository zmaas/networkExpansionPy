import networkExpansionPy.lib as ne
import networkExpansionPy.parseAtlas as pa
import networkExpansionPy.parseKegg as pk
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
from webweb import Web

# sets the expansion to run with the ATLAS database
atlas = True

# initializing network and pruning
# if atlas = True, this step will initalize with the atlas network database
metabolism = ne.GlobalMetabolicNetwork(atlas)
# if atlas = True, this step will skip all pruning and only convert to irreversible 
metabolism.init_pruning(atlas)

# Read in the seed compounds and parse them
cpds = pd.read_csv("../networkExpansionPy/assets/compounds/seeds_minimal.csv")
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())

# Run metabolic expansion
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)

# This designates which network reaction and compound IDs to associate with
if atlas == False:
    # Get dict of product/reactant pairs associated with each reaction
    rxn_pairs, _ = pa.get_rxn_pairs()  # FIXME may be redundant with Josh's
else:
    # Get dict of product/reactant pairs associated with each reaction
    rxn_pairs, _ = pk.get_rxn_pairs()  # FIXME may be redundant with Josh's 

print("----------------------------------------------------------------------")
print("Compounds: " + str(len(ne_cpds)))
print("Reactions: " + str(len(ne_rxns)))
print("----------------------------------------------------------------------")

# Get the number of new compounds at each step
rxns_count_list = [len(x) for x in ne_rxns_list]
rxns_count_differences = list(
    np.array(rxns_count_list[1:]) - np.array(rxns_count_list[:-1])
)
rxns_count_stepratios = list(
    np.array(rxns_count_differences) / np.array(rxns_count_list[:-1])
)
rxns_count_logquots = list(
    np.log2(np.array(rxns_count_list[1:]) / np.array(rxns_count_list[:-1]))
)
# Plot the growth over time
fig, ax = plt.subplots(nrows=2, ncols=2)
fig.suptitle("Per-Iteration Growth on Minimal Expansion Set with KEGG")
ax[0, 0].plot(rxns_count_list)
ax[0, 0].set_xlabel("Expansion Iteration")
ax[0, 0].set_ylabel("Total Reactions")
ax[0, 1].plot(rxns_count_differences)
ax[0, 1].set_xlabel("Expansion Iteration")
ax[0, 1].set_ylabel("Reactions Added")
ax[1, 0].plot(rxns_count_stepratios)
ax[1, 0].set_xlabel("Expansion Iteration")
ax[1, 0].set_ylabel("Reactions Added / Previous Network Size")
ax[1, 1].plot(rxns_count_logquots)
ax[1, 1].set_xlabel("Expansion Iteration")
ax[1, 1].set_ylabel("Log2(Relative Size Increase)")


# Read in data for name conversion
with open(
    "../networkExpansionPy/assets/iqbio/compounds.csv"
) as cpds_file_handle:
    lines = list(cpds_file_handle.readlines())
    cpd_dict = {
        cpd: name
        for cpd, name in [
            line.strip("cpd:").strip().split(maxsplit=1) for line in lines
        ]
    }

# Convert reactions to compounds
edges = set()
# Look up every product/reactant pair for all reactions in current step
for rxn in ne_rxns:
    try:
        edges.update(rxn_pairs[rxn[0]])
    except KeyError:
        print(f"Reaction {rxn[0]} not found")

try:      
    # Convert compounds to their real names (not CXXXXX)
    edges_renamed = [(cpd_dict[x], cpd_dict[y]) for x, y in edges]
    nodes_renamed = [cpd_dict[x] for x in ne_cpds]
    nodes_dct = {nodes_renamed[i]: i for i in range(0, len(nodes_renamed))}
    
    # Add a later to our webweb object
    web = Web(adjacency=list(edges_renamed), nodes=nodes_dct, title='Metabolism Network: Minimal')
    
    # Generate the webweb object for interactive visualization
    web.display.sizeBy = 'degree'
    web.display.colorBy = 'degree'
    web.show()
    
    # Code to identify seedset compounds that were not used in any reaction
    edges_set = set()
    for x,y in edges_renamed:
        edges_set.add(x)
        edges_set.add(y)
    
    list_of_unused = list()
    for k in nodes_renamed:
        if k not in edges_set:
            list_of_unused.append(k)
    print("Unused seedset compounds:")
    print(list_of_unused)

    
except KeyError:
    print("Compounds not found, cannot generate webweb")


# Code to create file for upload to KEGG
usr_data_kegg = open("../networkExpansionPy/assets/iqbio/Minimal_atlas.txt", "w")
for c in ne_cpds:
    usr_data_kegg.write(c + '\n')
# for r, direction in ne_rxns:
#     usr_data_kegg.write(r + '\n')

usr_data_kegg.close()

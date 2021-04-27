"""Runs a standard network expansion on our Europa seed set"""
import networkExpansionPy.lib as ne  # Expansion/graph code
import networkExpansionPy.parseKegg as pk  # Importing labels for plotting
import matplotlib.pyplot as plt  # For plotting
import numpy as np  # Numerical algorithms/data structures
import pandas as pd  # Data frames
import networkx as nx  # Graph algorithms/data structures
from webweb import Web  # Graph plotting

# initializing network and pruning
metabolism = ne.GlobalMetabolicNetwork()
metabolism.init_pruning()
metabolism.oxygen_independent()

# Read in the sead compounds and parse them
cpds = pd.read_csv("../networkExpansionPy/assets/compounds/seeds_europa.csv")
# Normalize by removing whitespace, using pandas formatting
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())

# Run metabolic expansion algorithm
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)

print("----------------------------------------------------------------------")
print("Compounds: " + str(len(ne_cpds)))
print("Reactions: " + str(len(ne_rxns)))

print("----------------------------------------------------------------------")

# Get the number of new compounds at each step
rxns_count_list = [len(x) for x in ne_rxns_list]
# Get the difference in compounds between subsequent steps
rxns_count_differences = list(
    np.array(rxns_count_list[1:]) - np.array(rxns_count_list[:-1])
)
# Calculate the relative network size increase
rxns_count_stepratios = list(
    np.array(rxns_count_differences) / np.array(rxns_count_list[:-1])
)
# Calculate the log2 ratio between subsequent steps
rxns_count_logquots = list(
    np.log2(np.array(rxns_count_list[1:]) / np.array(rxns_count_list[:-1]))
)

# # Plot the growth over time using the above data
# fig, ax = plt.subplots(nrows=2, ncols=2)
# fig.suptitle("Per-Iteration Growth on Europa Expansion Set")
# ax[0, 0].plot(rxns_count_list)
# ax[0, 0].set_xlabel("Expansion Iteration")
# ax[0, 0].set_ylabel("Total Reactions")
# ax[0, 1].plot(rxns_count_differences)
# ax[0, 1].set_xlabel("Expansion Iteration")
# ax[0, 1].set_ylabel("Reactions Added")
# ax[1, 0].plot(rxns_count_stepratios)
# ax[1, 0].set_xlabel("Expansion Iteration")
# ax[1, 0].set_ylabel("Reactions Added / Previous Network Size")
# ax[1, 1].plot(rxns_count_logquots)
# ax[1, 1].set_xlabel("Expansion Iteration")
# ax[1, 1].set_ylabel("Log2(Relative Size Increase)")


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

# Get dict of product/reactant pairs associated with each reaction
rxn_pairs, _ = pk.get_rxn_pairs()

# Convert reactions to compounds
edges = set()
# Look up every product/reactant pair for all reactions in current step
for rxn in ne_rxns:
    try:
        edges.update(rxn_pairs[rxn[0]])
    except KeyError:
        print(f"Reaction {rxn[0]} not found")

# Convert compounds to their real names (not CXXXXX)
edges_renamed = [(cpd_dict[x], cpd_dict[y]) for x, y in edges]
nodes_renamed = [cpd_dict[x] for x in ne_cpds]
nodes_dct = {nodes_renamed[i]: i for i in range(0, len(nodes_renamed))}

print("--------------------------------------------------")
print(ne_cpds)
print("--------------------------------------------------")


# Add a later to our webweb object
web = Web(adjacency=list(edges_renamed), nodes=nodes_dct, title='Metabolism Network: Europa')

# Generate the webweb object for interactive visualization
web.display.sizeBy = 'degree'
web.display.colorBy = 'degree'
web.show()


# Basic graph structure for exploring data
# G = nx.Graph()
# G.add_edges_from(edges_renamed)
# nx.drawing.nx_agraph.write_dot(G, "/hdd/minimal.dot")

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

# Code to create file for upload to KEGG
usr_data_kegg = open("../networkExpansionPy/assets/iqbio/Europa.txt", "w")
for c in ne_cpds:
    usr_data_kegg.write(c + '\n')
for r, direction in ne_rxns:
    usr_data_kegg.write(r + '\n')

usr_data_kegg.close()

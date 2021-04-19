import networkExpansionPy.lib as ne
import networkExpansionPy.parseKegg as pk
import networkExpansionPy.parseAtlas as pa
import pandas as pd
from webweb import Web
import matplotlib.pyplot as plt
import numpy as np

# initializing network and pruning
metabolism = ne.GlobalMetabolicNetwork(atlas=True)
metabolism.init_pruning(atlas=True)
metabolism.oxygen_indepentend()

# define seed compounds
cpds = pd.read_csv("../networkExpansionPy/assets/compounds/minimal_oxidized.csv")
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())

# run network expansion
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)

print("----------------------------------------------------------------------")
print("Compounds: " + str(len(ne_cpds)))
print("Reactions: " + str(len(ne_rxns)))

print("----------------------------------------------------------------------")
#print(ne_cpds)
#print("----------------------------------------------------------------------")
#print(ne_rxns)
#print("----------------------------------------------------------------------")

# plotting
rxn_pairs, _ = pa.get_rxn_pairs()

edges = set()

for rxn in ne_rxns:
    try:
        edges.update(rxn_pairs[rxn[0]])
    except KeyError:
        print(f"Reaction {rxn[0]} not found")

#print("----------------------------------------------------------------------")
#print(edges)

with open("../networkExpansionPy/assets/iqbio/compounds.csv") as cpds_file_handle:
    lines = list(cpds_file_handle.readlines())
    cpd_dict = {
        cpd: name
        for cpd, name in [
            line.strip("cpd:").strip().split(maxsplit=1) for line in lines
        ]
    }

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
fig.suptitle("Per-Iteration Growth on Minimal Expansion Set")
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

edges_renamed = [(cpd_dict[x], cpd_dict[y]) for x, y in edges]
nodes_renamed = [cpd_dict[x] for x in ne_cpds]
nodes_dct = {nodes_renamed[i]: i for i in range(0, len(nodes_renamed))}

#print("Nodes:")
#print(nodes_dct)
#print(type(nodes_dct))
#print("----------------------------------------------------------------------")
#print(edges_renamed)

web = Web(adjacency=list(edges_renamed), nodes=nodes_dct, title='Metabolism Network')
web.display.sizeBy = 'degree'
web.display.colorBy = 'degree'
web.show()

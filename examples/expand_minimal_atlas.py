import networkExpansionPy.lib as ne
import networkExpansionPy.parseAtlas as pa
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import networkx as nx
from webweb import Web

metabolism = ne.GlobalMetabolicNetwork(atlas=True)
# metabolism.pruneUnbalancedReactions()
# metabolism.pruneInconsistentReactions()
# metabolism.set_ph(7.0)
# metabolism.setMetaboliteBounds(ub=1e-1, lb=1e-6)
# metabolism.pruneThermodynamicallyInfeasibleReactions(keepnan=False)
metabolism.convertToIrreversible()

# Filter oxygen dependent reactions
oxygen_dependent_rxns = (
    metabolism.network[metabolism.network.cid.isin(["C00007"])]
    .rn.unique()
    .tolist()
)
o2_independent_rxns = [
    x
    for x in metabolism.network.rn.unique().tolist()
    if x not in oxygen_dependent_rxns
]
metabolism.subnetwork(o2_independent_rxns)

# Read in the sead compounds and parse them
cpds = pd.read_csv("../networkExpansionPy/assets/compounds/seeds_minimal.csv")
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())

# Run metabolic expansion
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)
# Get dict of product/reactant pairs associated with each reaction
rxn_pairs, _ = pa.get_rxn_pairs()  # FIXME may be redundant with Josh's

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

web = Web(title="Iterative Expansion from Minimal Seed Set")
for curr_rxns in ne_rxns_list:
    # Convert reactions to compounds
    edges = set()
    for rxn in curr_rxns:
        try:
            edges.update(rxn_pairs[rxn[0]])
        except KeyError:
            print(f"Lookup error on {rxn[0]}")
    # Convert compounds to their real names (not CXXXXX)
    edges_renamed = []
    for x, y in edges:
        try:
            x_name = cpd_dict[x]
        except KeyError:
            x_name = x
        try:
            y_name = cpd_dict[y]
        except KeyError:
            y_name = y
        edges_renamed.append((x_name,y_name))

    web.networks.webweb.add_layer(adjacency=edges_renamed)
# web.show()

# Basic graph structure for exploring data
# G = nx.Graph()
# G.add_edges_from(edges_renamed)
# nx.drawing.nx_agraph.write_dot(G, "/hdd/minimal.dot")

# Generate the webweb object for interactive visualization
# web.show()

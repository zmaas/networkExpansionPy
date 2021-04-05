"""Runs a standard network expansion on our Europa seed set"""
import networkExpansionPy.lib as ne  # Expansion/graph code
import networkExpansionPy.parseKegg as pk  # Importing labels for plotting
import matplotlib.pyplot as plt  # For plotting
import numpy as np  # Numerical algorithms/data structures
import pandas as pd  # Data frames
import networkx as nx  # Graph algorithms/data structures
from webweb import Web  # Graph plotting

# Initial metabolic network loading and pruning
# Construct network
metabolism = ne.GlobalMetabolicNetwork()
# Remove things with nonsense stoichiometries
metabolism.pruneUnbalancedReactions()
# Remove reactions that unrealistically produce new elements
metabolism.pruneInconsistentReactions()
# Look at a pH of 7 (must be between 5 and 9, 0.5 increments)
metabolism.set_ph(7.0)
# Upper and lower metabolite bounds
metabolism.setMetaboliteBounds(ub=1e-1, lb=1e-6)
# Irreversible required for thermodynamic considerations
metabolism.convertToIrreversible()
# Remove infeasible reactions
metabolism.pruneThermodynamicallyInfeasibleReactions(keepnan=False)

# Filter oxygen dependent reactions
# Shouldn't have an effect if O2 not present in seed?
filter_o2 = True
if filter_o2:
    # Find reactions requiring oxygen
    oxygen_dependent_rxns = (
        metabolism.network[metabolism.network.cid.isin(["C00007"])]
        .rn.unique()
        .tolist()
    )
    # Remove oxygen dependent reactions
    o2_independent_rxns = [
        x
        for x in metabolism.network.rn.unique().tolist()
        if x not in oxygen_dependent_rxns
    ]
    # Apply the subset to our network
    metabolism.subnetwork(o2_independent_rxns)

# Read in the sead compounds and parse them
cpds = pd.read_csv("../networkExpansionPy/assets/compounds/seeds_europa.csv")
# Normalize by removing whitespace, using pandas formatting
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())

# Run metabolic expansion algorithm
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)
# Get dict of product/reactant pairs associated with each reaction
rxn_pairs, _ = pk.get_rxn_pairs()  # FIXME may be redundant with Josh's

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
# Plot the growth over time using the above data
fig, ax = plt.subplots(nrows=2, ncols=2)
fig.suptitle("Per-Iteration Growth on Europa Expansion Set")
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

# Build webweb object
web = Web(title="Iterative Expansion from Minimal Seed Set")
# Grab the reactions and compound/compound pairs in each reaction
# FIXME the error we're getting can be fixed using the network df object, might be slower
for curr_rxns in ne_rxns_list:
    # Convert reactions to compounds
    edges = set()
    # Look up every product/reactant pair for all reactions in current step
    for rxn in curr_rxns:
        try:
            edges.update(rxn_pairs[rxn[0]])
        except KeyError:
            print(f"Reaction {rxn[0]} not found")
    # Convert compounds to their real names (not CXXXXX)
    edges_renamed = [(cpd_dict[x], cpd_dict[y]) for x, y in edges]
    # Add a later to our webweb object
    web.networks.webweb.add_layer(adjacency=edges_renamed)
# Generate the webweb object for interactive visualization
# web.show()

# Basic graph structure for exploring data
# G = nx.Graph()
# G.add_edges_from(edges_renamed)
# nx.drawing.nx_agraph.write_dot(G, "/hdd/minimal.dot")

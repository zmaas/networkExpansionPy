import networkExpansionPy.lib as ne
import networkExpansionPy.parseKegg as pk
import pandas as pd
import networkx as nx
from webweb import Web
import numpy as np

metabolism = ne.GlobalMetabolicNetwork()
metabolism.pruneUnbalancedReactions()
metabolism.pruneInconsistentReactions()
# metabolism.set_ph(7.0)
metabolism.convertToIrreversible()
# metabolism.setMetaboliteBounds(ub=1e-1,lb=1e-6)
# metabolism.pruneThermodynamicallyInfeasibleReactions(keepnan=False)

metabolism.convertToIrreversible()
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
# only keep anaerobic reactions
metabolism.subnetwork(o2_independent_rxns)
len(metabolism.network.cid.unique().tolist())
# define seed compounds
cpds = pd.read_csv("../networkExpansionPy/assets/compounds/seeds_minimal.csv")
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)
rxn_pairs, _ = pk.get_rxn_pairs()

edges = set()
for rxn in ne_rxns:
    try:
        edges.update(rxn_pairs[rxn[0]])
    except KeyError:
        print(f"Reaction {rxn[0]} not found")

with open("../networkExpansionPy/assets/iqbio/compounds.csv") as cpds_file_handle:
    lines = list(cpds_file_handle.readlines())
    cpd_dict = {
        cpd: name
        for cpd, name in [
            line.strip("cpd:").strip().split(maxsplit=1) for line in lines
        ]
    }

edges_renamed = [(cpd_dict[x], cpd_dict[y]) for x, y in edges]

# G = nx.Graph()
# G.add_edges_from(edges_renamed)
# nx.drawing.nx_agraph.write_dot(G, "../networkExpansionPy/assets/iqbio/compounds.csv")

web = Web(edges_renamed)
web.show()


print('Total compounds: {} Total reactions: {}'.format(np.size(ne_cpds), np.size(ne_rxns)/2))

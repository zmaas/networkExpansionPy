import networkExpansionPy.lib as ne
import networkExpansionPy.parseKegg as pk
import pandas as pd
from webweb import Web

# initializing network and pruning
metabolism = ne.GlobalMetabolicNetwork()
metabolism.init_pruning()
metabolism.oxygen_indepentend()

# define seed compounds
cpds = pd.read_csv("../networkExpansionPy/assets/iqbio/minimal_reduced.csv")
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
rxn_pairs, _ = pk.get_rxn_pairs()

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

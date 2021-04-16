import networkExpansionPy.lib as ne

# from scipy.sparse import csr_matrix
import numpy as np
# from random import sample
import pandas as pd

# from copy import copy, deepcopy
# import random
# import seaborn as sns
# from matplotlib_venn import venn2

df = pd.read_csv("../networkExpansionPy/assets/KEGG/rn2metal.csv")
rns = df[
    df["metal"].isin(
        ["fe", "fes", "fes_4Fe_4S", "fes_2Fe_2S", "fes_3Fe_4S", "mn"]
    )
].rxns.tolist()
metabolism = ne.GlobalMetabolicNetwork()
# metabolism.pruneUnbalancedReactions()
# metabolism.pruneInconsistentReactions()
# metabolism.set_ph(7.0)
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
metabolism.subnetwork(rns)

# metabolism.setMetaboliteBounds(ub=1e-1,lb=1e-6)
# metabolism.pruneThermodynamicallyInfeasibleReactions(keepnan=False)
len(metabolism.network.cid.unique().tolist())
# define seed compounds
cpds = pd.read_csv("../networkExpansionPy/assets/compounds/seeds.csv")
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)

print('Total compounds: {} Total reactions: {}'.format(np.size(ne_cpds), np.size(ne_rxns)/2))

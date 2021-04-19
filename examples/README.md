# Network Expansion Code Guide

<!-- Formatting: One sentence per line to make git diffing easier -->
## Running Network Expansion Code

TODO: UPDATE with current code

This is a quick guide on how to use the network expansion code in this repository to explore the expansion of your own seed set.
We'll explore how to construct a metabolic network, establish constraints, and then run a network expansion.
The first step is to load in the networkExpansionPy library:
```python
import networkExpansionPy.lib as ne  # Expansion/graph code
```
With the library loaded in we can construct the graph, and perform a few pruning steps:
```python
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
```
There are a few steps here that must be executed in order:
- Pruning inconsistent and unbalanced reactions should be first, to speed up execution time
- You must set the pH and metabolite bounds before pruning thermodynamically infeasible reactions
- You also must convert reactions to irreversible to remove infeasible reactions (delta G is direction dependent in the pruning process)

Once this is set up, you have a few options for how to proceed.
One option that is often useful is to remove reactions that depend on oxygen.
This can be done as follows:
```python
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
```

Whether or not you choose to remove oxygen dependent reactions, the next step is to load in your seed set.
This is a set of chemicals that you speculate are initially present that we will grow the network from.
The seed set is stored in a comma separated file (.csv) as follows:
```
CID,Name
C00001,Water
C00011,Carbon dioxide
C00080,Hydrogen Ions
C00014,Ammonia
C00009,Phosphate
C00283,Hydrogen sulfide
C01328,HO-
```
The most important thing in this file is making sure that the compound identifier (CID) is correct, as this is what the code will use internally when looking for reactions.
If this file was called `seeds.csv`, we could load it in as follows (using a full filesystem path if necessary:)
```python
# Read in the sead compounds and parse them
cpds = pd.read_csv("seeds.csv")
# Normalize by removing whitespace, using pandas formatting
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
seedset = set(cpds["CID"].tolist())
```
This code also does some formatting to remove any stray whitespace that might cause issues when looking up compound names.
Once we have this done, it's straightforward to run the actual network expansion:
```python
ne_cpds, ne_rxns, ne_cpds_list, ne_rxns_list = metabolism.expand(seedset)
```
This will return 4 variables:
- A list of compounds in the final expanded set
- A list of reactions in the final expanded set
- A list of lists of compounds at each step in the expansion
- A list of lists of reactions at each step in the expansion

Additionally a .txt file will be created of compound and reaction IDs.

These can then be used for later analysis and visualization.

## Visualization

TODO: add info on WebWeb and KEGG visualization steps with examples


## Functional Test

To check for functionality of the example expansions included, run the following test:
```
bash expansion_functional_tests.sh
```

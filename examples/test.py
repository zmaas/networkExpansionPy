# Import of the library file in networkExpansionPy
import networkExpansionPy.lib as ne

# Class scipy.sparse.csr_matrix
# (https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html#scipy.sparse.csr_matrix)
# A special format for matricies that mostly contaions 0s
# Sparse data are a data sets that mostly contains 0s
# CSR = Compressed Sparse Row
# Tutorial: https://www.w3schools.com/python/scipy_sparse_data.asp
from scipy.sparse import csr_matrix

# Python library to work with arrays
# https://www.w3schools.com/python/numpy_intro.asp
# In Python we have lists that serve the purpose of arrays, but they are slow to process
# NumPy aims to provide an array object that is up to 50x faster than traditional Python lists
import numpy as np

# Built-in module to generate "randomness"
# Tutorial for the sample function:
# https://www.w3schools.com/python/ref_random_sample.asp
# The sample() method returns a list with a randomly selection of a specified number of items from a sequnce
# random.sample(sequence, k)
from random import sample
import random

# Open source data analysis and manipulation tool
# https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html
import pandas as pd

# Data visualization library based on matplotlib (might not need it since
# Zach has ideas to visualize)
import seaborn as sns
from matplotlib_venn import venn2


# Reading comma seperated file into dataFrame
df = pd.read_csv('../networkExpansionPy/assets/KEGG/rn2metal.csv')

# Getting logical index of elements in metal column
index = df['metal'].isin(['fe','fes','fes_4Fe_4S','fes_2Fe_2S','fes_3Fe_4S','mn'])

# Converting corresponding rxns values into list
rns = df[index].rxns.tolist()

# Define metablism as class GlobalMetabolicNetwork
metabolism = ne.GlobalMetabolicNetwork()

# Modify metabolism
# metabolism.pruneUnbalancedReactions()
# metabolism.pruneInconsistentReactions()
# metabolism.set_ph(7.0)

# Split reversible reactions in two reactions (One forward and one backward)
metabolism.convertToIrreversible()

# Logical indexing of C00007 in column cid in network
# Storing the unique rn values in a list (rn is reaction identifier)
oxygen_dependent_rxns = metabolism.network[metabolism.network.cid.isin(['C00007'])].rn.unique().tolist()

# special python trick syntax to write small for loops that create a list
# add each unique element in rn to list if not part of O2 dependent list
o2_independent_rxns = [x for x in metabolism.network.rn.unique().tolist() if x not in oxygen_dependent_rxns]

# function that sets network to only contaion reactions that are in the list
# in this case only keep anaerobic reactions
metabolism.subnetwork(o2_independent_rxns)

# Reduse the network another time to only get reactions that correspond with the metals
metabolism.subnetwork(rns)

#metabolism.setMetaboliteBounds(ub=1e-1,lb=1e-6)
#metabolism.pruneThermodynamicallyInfeasibleReactions(keepnan=False)

# Get the number of unique cid values of the reduced network
# cid is compound ID
len(metabolism.network.cid.unique().tolist())

# At this point we shrunken the whole metabolic network to subnetwork that only includes
# reactions that are oxygen independent and reations that use the listed metals

# define seed compounds
# dataFrame that incldues CID and Name of the Seed components
cpds = pd.read_csv('/Users/joshuagoldford/Documents/github/networkExpansionPy/networkExpansionPy/assets/compounds/seeds.csv')

# A lambda function is a small anonymous function
# https://www.w3schools.com/python/python_lambda.asp
# lambda arguments : expression
# strip() removes the leading and trailing characters ???
cpds['CID'] = cpds['CID'].apply(lambda x: x.strip())

# Sets are used to store multiple items in a single variable
# https://www.w3schools.com/python/python_sets.asp
# Store the seedset
seedset = set(cpds['CID'].tolist())

# Run the expantion of the metabolism
# Returns the compounds and reactions
ne_cpds,ne_rxns = metabolism.expand(seedset)

# Give the number of compounds and reactions
len(ne_cpds)
len(seedset)

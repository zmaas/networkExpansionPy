##########################################################
#
# Team Kegger Network Expension Script
# Clair Huffine, Lukas Buecherl, Zachary Maas
#
# Advisors: Boswell Wing, Alexis S. Templeton
#
##########################################################

# Import of the library file in networkExpansionPy
import NetworkEuropa.lib as ne
# Open source data analysis and manipulation tool
# https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html
import pandas as pd
# Python library to work with arrays
# https://www.w3schools.com/python/numpy_intro.asp
# In Python we have lists that serve the purpose of arrays, but they are slow to process
# NumPy aims to provide an array object that is up to 50x faster than traditional Python lists
import numpy as np

# Define metablism as class GlobalMetabolicNetwork
metabolism = ne.GlobalMetabolicNetwork()

# Split reversible reactions in two reactions (One forward and one backward)
metabolism.convertToIrreversible()

# define seed compounds
# dataFrame that incldues CID and Name of the Seed components
seed = pd.read_csv('SeedEuropa.csv')

# Sets are used to store multiple items in a single variable
# https://www.w3schools.com/python/python_sets.asp
# Store the seedset
seedset = set(seed['CID'].tolist())

# Run the expantion of the metabolism
# Returns the compounds and reactions
ne_cpds,ne_rxns = metabolism.expand(seedset)

# Give the number of compounds and reactions
len(ne_cpds)
len(seedset)

import pandas as pd  # Data frames

europa = ['C00993', 'C00149', 'C00058', 'C00122', 'C00009', 'C00152', 'C01596', 'C00027', 'C00001', 'C19779', 'C14818', 'C00209', 'C00133', 'C01417', 'C01563', 'C00288', 'C00099', 'C00132', 'C00041', 'C00087', 'C00014', 'C00059', 'C00049', 'C00011', 'C00282', 'C00086', 'C00283', 'C00037', 'C00497', 'C01384', 'C00080', 'C00013', 'C00067', 'C00402']

# Read in the sead compounds and parse them
cpds = pd.read_csv("../networkExpansionPy/assets/iqbio/Target_Set.csv")
# Normalize by removing whitespace, using pandas formatting
cpds["CID"] = cpds["CID"].apply(lambda x: x.strip())
target = cpds["CID"].tolist()

set_europa = set(europa)
set_target = set(target)

print("-------------------------------")
print(set_target.intersection(set_europa))
print("-------------------------------")

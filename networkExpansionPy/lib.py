# Class scipy.sparse.csr_matrix
# (https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.csr_matrix.html#scipy.sparse.csr_matrix)
# A special format for matricies that mostly contaions 0s
# Sparse data are a data sets that mostly contains 0s
# CSR = Compressed Sparse Row
# Tutorial: https://www.w3schools.com/python/scipy_sparse_data.asp
from scipy.sparse import csr_matrix
from scipy.sparse import dok_matrix

# Python library to work with arrays
# https://www.w3schools.com/python/numpy_intro.asp
# In Python we have lists that serve the purpose of arrays, but they are slow to process
# NumPy aims to provide an array object that is up to 50x faster than traditional Python lists
import numpy as np

# Open source data analysis and manipulation tool
# https://pandas.pydata.org/pandas-docs/stable/user_guide/10min.html
import pandas as pd

# Python package that helps interating with your operating system
import os

# Python package to work with JSON data. JSON is a syntax for storing and exchanging data, it is
# text, written with JavaScript object notation.
import json

# Package that adds copy functions to python. Python usually uses references, but
# sometimes real copies are necessary to edit one copy without changing the other
from copy import copy, deepcopy

# Progress bars
from tqdm import tqdm

# Define the asset path
# The first line gets the path of this file itself and splits the path into head and tail.
# Tail is the last named component (here: lib.py) head is the rest.
asset_path,filename = os.path.split(os.path.abspath(__file__))
# The asset path is defined by using the head of the previous function and adding /asset
asset_path = asset_path + '/assets'

# Algorithm that runs the network expansion. More information can be found in the method
# section of Josh's paper or in our notes document: https://docs.google.com/document/d/17I7lOJJPZ8XeVjlQ45GArG3s_Yey-CnGIgjZFpHG0Os/edit
def netExp(R, P, x, b):
    k = np.sum(x)
    k0 = 0
    n_reactions = np.size(R, 1)
    y = csr_matrix((n_reactions))
    x_list = [x]
    y_list = [y]
    pbar = tqdm(desc="Network Expansion")
    while k > k0:
        k0 = np.sum(x)
        y = np.dot(R.transpose(), x) == b
        y = y.astype("int")
        x_n = np.dot(P, y) + x
        x_n = x_n.astype("bool")
        x = x_n.astype("int")
        k = np.sum(x)
        x_list.append(x)
        y_list.append(y)
        pbar.update(1)
    pbar.close()
    return x, y, x_list, y_list

# Algorithm that runs the network expansion, here the stopping criteria is not only that no new compounds are added,
# but also that no new reactions are added. (Another counter l is added to keep track of reactions)
def netExp_cr(R,P,x,b):
    k = np.sum(x);
    k0 = 0;
    n_reactions = np.size(R,1)
    y = csr_matrix(np.zeros(n_reactions))
    l = 0
    l0 = 0;

    while (k > k0) | (l > l0):
        k0 = np.sum(x);
        l0 = np.sum(y)
        y = (np.dot(R.transpose(),x) == b);
        y = y.astype('int');
        x_n = np.dot(P,y) + x;
        x_n = x_n.astype('bool');
        x = x_n.astype('int');
        k = np.sum(x);
        l = np.sum(y)
    return x,y


# Algorithm that runs the network expansion with no new compounds as stopping criteria
# Here, lists are defined and after each iteration, x and y is added to the list
# List items are ordered, changeable, and allow duplicate values.
def netExp_trace(R,P,x,b):

    # Define lists
    X = []
    Y = []

    X.append(x)
    k = np.sum(x);
    k0 = 0;
    n_reactions = np.size(R,1)
    y = csr_matrix(np.zeros(n_reactions))
    Y.append(y)

    while k > k0:
        k0 = np.sum(x);
        y = (np.dot(R.transpose(),x) == b);
        y = y.astype('int');
        x_n = np.dot(P,y) + x;
        x_n = x_n.astype('bool');
        x = x_n.astype('int');
        k = np.sum(x);
        X.append(x)
        Y.append(y)
    return X,Y

# ???
def parse_reaction_trace(reaction_trace,network):
    rxns_list = []
    for i in range(1,len(reaction_trace)):
        idx = reaction_trace[i].nonzero()[0]
        rxns = list(network.iloc[:,idx])
        rxns = pd.DataFrame(rxns,columns = ['rn','direction'])
        rxns['iter'] = i
        rxns_list.append(rxns)
    rxns_list = pd.concat(rxns_list,axis=0)
    return rxns_list


# ???
def isRxnCoenzymeCoupled(rxn,cosubstrate,coproduct):
    g = rxn[rxn.cid.isin([cosubstrate,coproduct])]
    out = False
    if len(g) > 1:
        if g.s.sum() == 0:
            out = True
    return out


# ???
def load_ecg_network(ecg):
    network_list = []
    consistent_rids = []
    for rid,v in ecg["reactions"].items():
        cids = v["left"] + v["right"]
        try: ## This skips all reactions with n stoichiometries
            stoichs = [-int(i) for i in v["metadata"]["left_stoichiometries"]]+[int(i) for i in v["metadata"]["right_stoichiometries"]]
            network_list+=list(zip(cids,[rid for _ in range(len(stoichs))],stoichs))
        except:
            pass
        if v["metadata"]["element_conservation"]==True:
            consistent_rids.append(rid)
    return pd.DataFrame(network_list,columns=("cid","rn","s")), pd.DataFrame(consistent_rids,columns=["rn"])


# The class the defines the characteristics of the metabolic network
class GlobalMetabolicNetwork:

    # The init function gets called when a class object is initialized
    # the function checks whether or not to load the atlas or kegg database with 
    # the kegg database as default.
    # The function checks if a databse in json format is given, if not it loads
    # the /KEGG/network_full.csv' as default. It also loads /compounds/cpds.txt'
    # for the compounds and /reaction_free_energy/kegg_reactions_CC_ph7.0.csv'
    # for the thermodynamics
    def __init__(self,atlas=None, ecg_json=None):

        if ecg_json == None:
            if atlas == None:
                network = pd.read_csv(asset_path + '/KEGG/network_full.csv')
            else:
                network = pd.read_csv(asset_path + '/KEGG/atlas_network_full.csv')
                print('Running with ATLAS database')
            cpds = pd.read_csv(asset_path +'/compounds/cpds.txt',sep='\t')
            thermo = pd.read_csv(asset_path +'/reaction_free_energy/kegg_reactions_CC_ph7.0.csv',sep=',')
            self.network = network
            self.thermo = thermo
            self.compounds = cpds ## Includes many compounds without reactions
            self.ecg = None
        else:
        # If other network is specified use the other network instead
            with open(ecg_json) as f:
                ecg = json.load(f)
            network, consistent_rxns = load_ecg_network(ecg)
            self.network = network
            self.ecg = ecg
            self.consistent_rxns = consistent_rxns
            self.compounds = pd.DataFrame(self.network["cid"].unique(),columns=["cid"]) ## Only iitncludes compounds with reactions

        # Set attributes of the class (Temp is initialized to 25C)
        self.temperature = 25
        self.seedSet = None
        self.rid_to_idx = None
        self.idx_to_rid = None
        self.cid_to_idx = None
        self.idx_to_cid = None
        self.S = None

    # Function that uses the copy library and returns a true copy of the class object
    # (Not just a reference)
    def copy(self):
        return deepcopy(self)

    # Function to set the pH of the network. If no file is specifeid it uses '/reaction_free_energy/kegg_reactions_CC_ph'
    # by default
    def set_ph(self,pH):
        if ~(type(pH) == str):
            pH = str(pH)
        if self.ecg == None:
            try:
                thermo = pd.read_csv(asset_path + '/reaction_free_energy/kegg_reactions_CC_ph' + pH + '.csv',sep=',')
                self.thermo = thermo
            except Exception as error:
                print('Failed to open pH files (please use 5.0-9.0 in 0.5 increments)')
        else:
            try:
                self.thermo = load_ecg_thermo(self.ecg,pH)
            except:
                raise ValueError("Try another pH, that one appears not to be in the ecg json")

    # ???
    def load_ecg_thermo(self,ph=9):
        thermo_list = []
        for rid,v in self.ecg["reactions"].items():

            phkey = str(ph)+"pH_100mM"

            if v["metadata"]["dg"][phkey]["standard_dg_prime_value"] == None:
                dg = np.nan
            else:
                dg = v["metadata"]["dg"][phkey]["standard_dg_prime_value"]

            if v["metadata"]["dg"][phkey]["standard_dg_prime_error"] == None:
                dgerror = np.nan
            else:
                dgerror = v["metadata"]["dg"][phkey]["standard_dg_prime_error"]

            if v["metadata"]["dg"][phkey]["is_uncertain"] == None:
                note = "uncertainty is too high"
            else:
                note = np.nan

            thermo_list.append((rid,
                dg,
                dgerror,
                v["metadata"]["dg"][phkey]["p_h"],
                v["metadata"]["dg"][phkey]["ionic_strength"]/1000,
                v["metadata"]["dg"][phkey]["temperature"],
                note))

        return pd.DataFrame(thermo_list, columns = ("!MiriamID::urn:miriam:kegg.reaction","!dG0_prime (kJ/mol)","!sigma[dG0] (kJ/mol)","!pH","!I (mM)","!T (Kelvin)","!Note"))

    # This function removes reactions that produce new elements by checking if the elements in
    # the reaction and products are the same
    def pruneInconsistentReactions(self):
        # remove reactions with qualitatively different sets of elements in reactions and products
        if self.ecg==None:
            consistent = pd.read_csv(asset_path + '/reaction_sets/reactions_consistent.csv')
            self.network = self.network[self.network.rn.isin(consistent.rn.tolist())]
        else:
            self.network = self.network[self.network.rn.isin(self.consistent_rxns.rn.tolist())]

    # This funciton removes reactions with inconsistent stoichiometries by comparing
    # the network to '/reaction_sets/reactions_balanced.csv'
    def pruneUnbalancedReactions(self):
        # only keep reactions that are elementally balanced
        balanced = pd.read_csv(asset_path + '/reaction_sets/reactions_balanced.csv')
        self.network = self.network[self.network.rn.isin(balanced.rn.tolist())]

    # This function creates a subnetwork of the given network by comparing and only keeping
    # elements that are specified in the given list
    def subnetwork(self,rxns):
        # only keep reactions that are in list
        self.network = self.network[self.network.rn.isin(rxns)]

    # This function replaces NAD/NADP/FAD with generic coenzymes
    def addGenericCoenzymes(self):
        replace_metabolites = {'C00003': 'Generic_oxidant', 'C00004': 'Generic_reductant', 'C00006': 'Generic_oxidant',  'C00005': 'Generic_reductant','C00016': 'Generic_oxidant','C01352':'Generic_reductant'}
        coenzyme_pairs = {}
        coenzyme_pairs['NAD'] = ['C00003','C00004']
        coenzyme_pairs['NADP'] = ['C00006','C00005']
        coenzyme_pairs['FAD'] = ['C00016','C01352']
        coenzyme_pairs = pd.DataFrame(coenzyme_pairs).T.reset_index()
        coenzyme_pairs.columns = ['id','oxidant','reductant']
        # create reactions copies with coenzyme pairs
        new_rxns = []
        new_thermo = [];
        for idx,rxn in self.network.groupby('rn'):
            z = any([isRxnCoenzymeCoupled(rxn,row.oxidant,row.reductant) for x,row in coenzyme_pairs.iterrows()])
            if z:
                new_rxn = rxn.replace(replace_metabolites).groupby(['cid','rn']).sum().reset_index()
                new_rxn['rn'] = new_rxn['rn'] = idx + '_G'
                new_rxns.append(new_rxn)
                t = self.thermo[self.thermo['!MiriamID::urn:miriam:kegg.reaction'] == idx].replace({idx:  idx + '_G'})
                new_thermo.append(t)

        new_rxns = pd.concat(new_rxns,axis=0)
        new_thermo = pd.concat(new_thermo,axis=0)

        self.network = pd.concat([self.network,new_rxns],axis=0)
        self.thermo = pd.concat([self.thermo,new_thermo],axis=0)


    # This function converts reversible reactions in irreversible reactions by
    # copying the network and multiplying one copy by -1 and labeling the copies
    # forward and reverse
    def convertToIrreversible(self):
        # Split reversible reactions into one forward and one backward reaction
        # copy the network
        nf = self.network.copy()
        nb = self.network.copy()

        nf['direction'] = 'forward'
        nb['direction'] = 'reverse'
        nb['s'] = -nb['s']

        # Concatenate nf,nb
        net = pd.concat([nf,nb],axis=0)
        net = net.set_index(['cid','rn','direction']).reset_index()
        self.network = net

    # This function sets the upper and lower bounds of the metabolits. By default
    # the upper bound is set to 1e-1 und the lower bound to 1e-6
    def setMetaboliteBounds(self,ub = 1e-1,lb = 1e-6):

        self.network['ub'] = ub
        self.network['lb'] = lb

    # This function calcualtes the effective delta G of each reaction and drops reactions
    # that are thermodynamical infeasable
    def pruneThermodynamicallyInfeasibleReactions(self,keepnan = False):
        fixed_mets = ['C00001','C00080']

        RT = 0.008309424 * (273.15+self.temperature)
        rns  = []
        dirs = []
        dgs = []
        for (rn,direction), dff in self.network.groupby(['rn','direction']):
            effective_deltaG = np.nan
            if rn in self.thermo['!MiriamID::urn:miriam:kegg.reaction'].tolist():
                deltaG = self.thermo[self.thermo['!MiriamID::urn:miriam:kegg.reaction'] == rn]['!dG0_prime (kJ/mol)'].values[0]
                if direction == 'reverse':
                    deltaG = -1*deltaG

                dff = dff[~dff['cid'].isin(fixed_mets)]
                subs = dff[dff['s'] < 0]
                prods = dff[dff['s'] > 0];
                k = np.dot(subs['ub'].apply(np.log),subs['s']) + np.dot(prods['lb'].apply(np.log),prods['s'])

                effective_deltaG = RT*k + deltaG

            dgs.append(effective_deltaG)
            dirs.append(direction)
            rns.append(rn)

        res = pd.DataFrame({'rn':rns,'direction':dirs,'effDeltaG':dgs})
        if ~keepnan:
            res = res.dropna()

        #res = res[res['effDeltaG'] < 0].set_index(['rn','direction'])
        res = res[~(res['effDeltaG'] > 0)].set_index(['rn','direction'])
        res = res.drop('effDeltaG',axis=1)
        self.network = res.join(self.network.set_index(['rn','direction'])).reset_index()


    # This function initializes the vector x. If there are n metabolits x has the
    # length of n. If reactant i is present xi is one otherwise it is 0
    def initialize_metabolite_vector(self,seedSet):
        if seedSet is None:
            print('No seed set')
        else:
            x0 = np.zeros([len(self.cid_to_idx)],dtype=int)
            for x in set(seedSet)&set(self.cid_to_idx.keys()):
                x0[self.cid_to_idx[x]] = 1
            return x0

    # This function creates dictionaries to store reactions. Dictionaries are used
    # to store data values in key:value pairs. Here doubles of (rn,direction) are
    # created and put into a set. Those sets are added to the dictionaries with
    # the key idx
    def create_reaction_dicts(self):
        # Create doubles of (rn,direction) and put them in a set
        rids = set(zip(self.network["rn"],self.network["direction"]))
        # Initialize dictionaries
        rid_to_idx = dict()
        idx_to_rid = dict()
        for v, k in enumerate(rids):
            rid_to_idx[k] = v
            idx_to_rid[v] = k

        return rid_to_idx, idx_to_rid

    # This function works like create_reaction_dicts(self). It creates dictionaries
    # to store compounds with a key (idx)
    def create_compound_dicts(self):
        cids = set(self.network["cid"])
        # Initialize dictionaries
        cid_to_idx = dict()
        idx_to_cid = dict()
        for v, k in enumerate(cids):
            cid_to_idx[k] = v
            idx_to_cid[v] = k

        return cid_to_idx, idx_to_cid

    # This function creates a matrix S that contains the compound id, reaction,
    # direction, and stoichiometrie of a metabolite
    def create_S_from_irreversible_network(self):

       S = dok_matrix((len(self.cid_to_idx),len(self.rid_to_idx)))

       for c,r,d,s in tqdm(zip(self.network["cid"],self.network["rn"],self.network["direction"],self.network["s"]), desc="Converting to Irreversible", total=len(self.network)):
           S[self.cid_to_idx[c],self.rid_to_idx[(r,d)]] = s

       return S.tocsr()

    # The network expanation algorithm that utilizes the netExp() function at the beginng of the
    # library. The function's inputs are the seedSet and the algorithm specification, specifying
    # if the normal or reaction and compound restricted algorithm should be ran. The function's
    # output are two lists one containing the final compounds and one containing the final reactions
    # of the expanded network.
    def expand(self,seedSet,algorithm='naive'):
        # constructre network from skinny table and create matricies for NE algorithm
        # if (self.rid_to_idx is None) or (self.idx_to_rid is None):
        self.rid_to_idx, self.idx_to_rid = self.create_reaction_dicts()
        # if (self.cid_to_idx is None) or (self.idx_to_cid is None):
        self.cid_to_idx, self.idx_to_cid = self.create_compound_dicts()
        # if self.S is None:
        self.S = self.create_S_from_irreversible_network()
        x0 = self.initialize_metabolite_vector(seedSet)
        R = (self.S < 0)*1
        P = (self.S > 0)*1
        b = sum(R)

        # sparsefy data
        R = csr_matrix(R)
        P = csr_matrix(P)
        b = csr_matrix(b)
        b = b.transpose()

        x0 = csr_matrix(x0)
        x0 = x0.transpose()
        if algorithm.lower() == "naive":
            x, y, x_list, y_list = netExp(R, P, x0, b)
        elif algorithm.lower() == "cr":
            x, y = netExp_cr(R, P, x0, b)
        else:
            raise ValueError(
                "algorithm needs to be naive (compound stopping criteria) or cr (reaction/compound stopping criteria)"
            )

        # convert to list of metabolite ids and reaction ids
        if x.toarray().sum() > 0:
            cidx = np.nonzero(x.toarray().T[0])[0]
            compounds = [self.idx_to_cid[i] for i in cidx]
        else:
            compounds = []

        if y.toarray().sum() > 0:
            ridx = np.nonzero(y.toarray().T[0])[0]
            reactions = [self.idx_to_rid[i] for i in ridx]
        else:
            reactions = []

        compounds_list = []
        for sub_x in x_list:
            if sub_x.toarray().sum() > 0:
                cid_sub_x = np.nonzero(sub_x.toarray().T[0])[0]
                sub_compounds = [self.idx_to_cid[i] for i in cid_sub_x]
            else:
                sub_compounds = []
            compounds_list.append(sub_compounds)

        reactions_list = []
        for sub_x in x_list:
            if sub_x.toarray().sum() > 0:
                cid_sub_x = np.nonzero(sub_x.toarray().T[0])[0]
                sub_reactions = [self.idx_to_rid[i] for i in cid_sub_x]
            else:
                sub_reactions = []
            reactions_list.append(sub_reactions)

        return compounds, reactions, compounds_list, reactions_list

    # Function to prune the network initially  
    def init_pruning(self,pH='7.0', ub=1e-1, lb=1e-6, keepnan=False):
        self.pruneUnbalancedReactions()
        # Remove reactions that unrealistically produce new elements
        self.pruneInconsistentReactions()
        # Look at a pH of 7 (must be between 5 and 9, 0.5 increments)
        self.set_ph(pH)
        # Upper and lower metabolite bounds
        self.setMetaboliteBounds(ub, lb)
        # Irreversible required for thermodynamic considerations
        self.convertToIrreversible()
        # Remove infeasible reactions
        self.pruneThermodynamicallyInfeasibleReactions(keepnan)

    # Function to delete all oxygen dependent reactions
    def oxygen_indepentend(self):
        oxygen_dependent_rxns = (
            self.network[self.network.cid.isin(["C00007"])]
            .rn.unique()
            .tolist()
        )
        o2_independent_rxns = [
            x
            for x in self.network.rn.unique().tolist()
            if x not in oxygen_dependent_rxns
        ]
        # only keep anaerobic reactions
        self.subnetwork(o2_independent_rxns)

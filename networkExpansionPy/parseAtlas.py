# parseKegg.py --- Parse downloaded KEGG data from getKegg.py
#
# Filename: parseAtlas.py
# Author: Clair Huffine <clair.huffine@colorado.edu>
# Created: Thu Apr 15 09:56:26 2021 (-0600)
#
#

# Commentary:
#
#
# This file implements a rudimentary parser to convert downloaded ATLAS
# files into a basic graph implementation. Here, we don't care about
# stoichiometries yet, and edges are defined as all product/reactant
# pairs. This isn't the best way to do this but it's easy.
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
#
#

# Code:

import re
import itertools as it
import networkx as nx
import matplotlib.pyplot as plt
import sys
from array import array
import csv


in_file = "../networkExpansionPy/assets/iqbio/atlas_full.csv"
out_file = "../networkExpansionPy/assets/iqbio/atlas_network_full.csv"

def get_rxn_pairs(in_file= "../networkExpansionPy/assets/iqbio/atlas_full.csv", init= 1 ):
    all_rxns = get_column(in_file, [0, 2])
    identify_compounds(all_rxns, init)

def get_column(file_name, result_column=1):

    """ Opens a comma separated CSV file and returns an array of intergers
        from the result_column,

    Parameters
    ----------
    file_name: string
            The path to the CSV file.
    result_column: integer or string or tuple
            The column to be returned.
            If a single column, works only on a string that can be made into
            an integer. If multiple columns, will work with ints only.

    Returns:
    --------
    Output: array of intergers or list of lists with multiple result column
            inputs values in the result_column.
    """

    # Checks for file not found and perrmission errors
    try:
        f = open(file_name, 'r')
    except FileNotFoundError:
        print("Couldn't find file " + file_name)
        sys.exit(3)
    except PermissionError:
        print("Couldn't access file " + file_name)
        sys.exit(4)

    f = open(file_name, 'r', encoding="ISO-8859-1")

    # separates header line, removes /n, and splits into its elements
    header = f.readline().rstrip().split(',')

    # calls the column_index function to identify the query and result columns
    # based on either their integer or string value
    column_index = identify_column(result_column, header)

    # result column index
    ii = column_index

    # sets Result output type based on number of inputted result columns
    if type(result_column) is list:
        # create Results list to add results to with multiple result columns
        Output = []
    else:
        Output = array('i', [])

    for line in f:
        A = line.rstrip().split(',')

        # appends value in the result columns to the outputted Result list
        if type(result_column) is list:
            case = []
            for column in result_column:
                case.append(A[int(column)])
            Output.append(case)
        # appends value in the result column to the outputted Result array
        else:
            Output.append(int(A[ii]))

    f.close()
    return(Output)

def identify_column(result_column, header):

    """ Identifies the result columns of interest if entered as
        an integer or a string and gives the index values of those columns.

    Parameters
    ----------
    result_column: integer or string or tuple
            The column to be returned.
            If a single column, works only on a string that can be made into
            an integer. If multiple columns, will work with ints only.
    header: list of strings
            The header of the file of interest which contains the column
            titles to be searched for a match with inputted strings for result
            columns.

    Returns:
    --------
    ii: integer
            Index value of the result column.
    """

    # generates variable ii to indicate the result column
    ii = 0
    if type(result_column) is not list:
        # checks if integer for column was inputted
        # assumes counting starting at 1
        for column_header2 in header:
            try:
                type(int(result_column)) == int
                result_column = int(result_column)
                # checks for array length exception
                if result_column > (len(header) - 1):
                    print("Searching for result column "
                          + str(result_column)
                          + " but there are only "
                          + str(len(header)-1)
                          + " fields")
                    sys.exit(1)
                else:
                    ii = result_column - 1
            except ValueError:
                if result_column == column_header2:
                    break
                else:
                    ii += 1
                    # checks for str(query_column) match exception
                    if ii > (len(header)-1):
                        print("Searching for result column "
                              + result_column
                              + " but this file does not contain it")
                        sys.exit(1)
    else:
        for column in result_column:
            # checks for array length exception
            if int(column) > (len(header)-1):
                print("Searching for result column "
                      + str(column)
                      + " but there are only "
                      + str(len(header)-1)
                      + " fields")
                sys.exit(1)
    return(ii)

def identify_compounds(all_rxns, init=0):
    
    """ Takes reactions ID and reaction equations and splits the reactants and
    products for each reaction. Outputs a dictonary of the ATLAS reaction ID 
    and a list of the reactants and products. 

    Parameters
    ----------
    all_rxns: list of lists of strings
            A list of lists containing each reaction ID and the associated 
            reaction from the ATLAS data set. 

    Returns:
    --------
    rxn_pairs: dictonary
            A dictonary of the ATLAS reaction ID and a list of the paired
            reactants and products.
    all_pairs: set
            Reactant/product pairs from each reaction
    """
    
    all_pairs = set()
    regex = r"C\d{5}"
    rxn_pairs = {}

    for rxn, eqn in all_rxns:
        # Break out and find reactions and products
        reactant_str, product_str = eqn.split("<=>")
        if init == 0:
            compounds = {}
            compounds[-1] = re.findall(regex, reactant_str)
            compounds[1] = re.findall(regex, product_str)
            pairs = list(it.product(compounds[-1], compounds[1]))
            rxn_pairs[rxn] = (compounds)
        elif init == 1:
            reactants = re.findall(regex, reactant_str)
            products = re.findall(regex, product_str)
            # Generate pairwise combinations as edges
            pairs = list(it.product(reactants, products))
            rxn_pairs[rxn] = pairs

        all_pairs.update(pairs)

    return (rxn_pairs, all_pairs)



def create_full_network(rxn_pairs):
       
    """ Takes a dictonary of the ATLAS reaction ID and a list of the reactants
    and products and writes them into a .csv file with each row consisting of 
    a compound ID (cid) and a reaction ID (rn) it is associated with. 

    Parameters
    ----------
    rxn_pairs: dictonary
            A dictonary of the ATLAS reaction ID and a list of the paired
            reactants and products.

    Outputs:
    --------
    out_file: .csv
            A .csv file with each row consisting of a compound ID (cid) and a
            reaction ID (rn) it is associated with.
    """

    with open(out_file, "w+", newline='') as rxn_file_handle:
        network_writer = csv.writer(rxn_file_handle, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        network_writer.writerow(('cid', 'rn', 's'))
        # Iterate over all reactions, writing result to file
        for rxn_ID in rxn_pairs:
            for reactant_ID in rxn_pairs[rxn_ID][-1]:
                network_writer.writerow((reactant_ID, rxn_ID, -1))
            for product_ID in rxn_pairs[rxn_ID][1]:
                network_writer.writerow((product_ID, rxn_ID, 1))

def main():

    all_rxns = get_column(in_file, [0, 2])
    rxn_pairs, all_pairs = identify_compounds(all_rxns)
    create_full_network(rxn_pairs)

if __name__ == '__main__':
    main()

#
# parseAtlas.py ends here

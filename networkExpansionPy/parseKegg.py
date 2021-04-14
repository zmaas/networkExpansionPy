# parseKegg.py --- Parse downloaded KEGG data from getKegg.py
#
# Filename: parseKegg.py
# Author: Zachary Maas <zama8258@colorado.edu>
# Created: Thu Mar 18 12:04:26 2021 (-0600)
#
#

# Commentary:
#
#
# This file implements a rudimentary parser to convert downloaded KEGG
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

in_file = "../networkExpansionPy/assets/iqbio/kegg_rxns.tab"
out_file = "../networkExpansionPy/assets/iqbio/kegg_rxns.dot"
out_gmlfile = "../networkExpansionPy/assets/iqbio/kegg_rxns.gml"


def get_rxn_pairs(in_file=in_file):
    all_pairs = set()
    with open(in_file, "r") as rxn_file_handle:
        lines = list(rxn_file_handle.readlines())
        # Grab only the equations from this list
        rxns = [
            line.strip("ENTRY").strip().split()[0]
            for line in lines
            if "ENTRY" in line
        ]
        eqns = [
            line.strip("EQUATION").strip()
            for line in lines
            if "EQUATION" in line
        ]
        if len(rxns) != len(eqns):
            raise Exception(
                "Unequal number of reactions and equations in data file."
            )
        regex = r"C\d{5}"
        rxn_pairs = {}
        for rxn, eqn in zip(rxns, eqns):
            # Break out and find reactions and products
            reactant_str, product_str = eqn.split("<=>")
            reactants = re.findall(regex, reactant_str)
            products = re.findall(regex, product_str)
            # Generate pairwise combinations as edges
            pairs = list(it.product(reactants, products))
            rxn_pairs[rxn] = pairs
            all_pairs.update(pairs)
    return (rxn_pairs, all_pairs)


def gen_filtered_graph(all_pairs):
    # Generate Graph
    G = nx.Graph()
    G.add_edges_from(all_pairs)

    # Filter by degree > 5 (arbitrary)
    deg = G.degree()
    to_keep = [n for n, d in deg if d > 20]
    sG = G.subgraph(to_keep)

    nx.drawing.nx_agraph.write_dot(sG, out_file)
    nx.write_graphml(G, out_gmlfile)


#
# parseKegg.py ends here

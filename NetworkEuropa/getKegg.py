# getKegg.py --- Get KEGG Database
#
# Filename: getKegg.py
# Author: Zachary Maas <zama8258@colorado.edu>
# Created: Thu Mar 18 11:31:57 2021 (-0600)
#
#

# Commentary:
#
#
# This script will pull the full KEGG metabolism pathway from the
# online API without having to use the FTP server. It would obviously
# be easier to use the FTP server but we're not rich so the API will
# have to do for now.
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

from itertools import zip_longest
import requests
from tqdm import tqdm

# Output file path, you should change this
out_file = "/hdd/kegg_rxns.tab"


def grouper(iterable, n, fillvalue=None):
    """Group ITERABLE into lists of size N, padding with fillvalue if
    len(ITERABLE) mod n is nonzero. Essentially, split any iterable
    (e.g. a list) into a list of lists of size n, padding the last
    list of lists with fillvalue if the size of iterable is not
    exactly divisible by n.

    """
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue=fillvalue)


# Session use gives us a speedup and reuses the same connection
with requests.Session() as http_sess:
    # Get the reaction list and parse out reaction identifiers
    rxns_response = http_sess.get("http://rest.kegg.jp/list/reaction")
    rxns = [
        rxn_str.split("\t")[0] for rxn_str in rxns_response.text.split("\n")
    ]
    # Chunk into groups of 10 per KEGG's limitations
    rxns_chunked = list(grouper(rxns, 10, ""))
    # Output file is user specified
    with open(out_file, "w+") as rxn_file_handle:
        # Iterate over all reactions, writing result to file
        for rxn in tqdm(rxns_chunked, desc="Network Requests"):
            rxn_response = http_sess.get(
                f"http://rest.kegg.jp/get/{'+'.join(rxn)}"
            )
            rxn_file_handle.write(rxn_response.text)


#
# getKegg.py ends here

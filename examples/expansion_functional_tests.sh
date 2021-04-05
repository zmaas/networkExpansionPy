
!/bin/bash

"""
Created on Mon Apr  5 14:09:40 2021

@author: Clair Huffine

This is a file of functional tests for network expansion to assure correct runs
across different expasions. To be run after alterations of foundational code in 
the lib.py file or elsewhere. 

This test relies on Ryan Layer's stupid simple ba(sh) testing found at:
https://github.com/ryanlayer/ssshtest
"""

test -e ssshtest || wget -q
https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run test_full_run expand_europa.py
assert_exit_code 0
# assert_in_stdout 
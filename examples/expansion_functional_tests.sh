#!/bin/bash

# @author: Clair Huffine
# This is a file of functional tests for network expansion to assure correct runs
# across different expasions. To be run after alterations of foundational code in 
# the lib.py file or elsewhere. 
# This test relies on Ryan Layer's stupid simple ba(sh) testing found at:
# https://github.com/ryanlayer/ssshtest

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run test_example python -c "print ('example assert_in_stdout success')"
assert_in_stdout "example"

run test_minimal python expand_minimal.py
assert_exit_code 0
assert_in_stdout "Total compounds: 29 Total reactions: 54"


run test_metals python expand_metals.py
assert_exit_code 0
assert_in_stdout "Total compounds: 56 Total reactions: 68"

run test_europa python expand_europa.py
assert_exit_code 0
assert_in_stdout "Total compounds: 32 Total reactions: 34"
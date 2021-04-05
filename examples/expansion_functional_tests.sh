#!/bin/bash

# @author: Clair Huffine
# This is a file of functional tests for network expansion to assure correct runs
# across different expasions. To be run after alterations of foundational code in 
# the lib.py file or elsewhere. 
# This test relies on Ryan Layer's stupid simple ba(sh) testing found at:
# https://github.com/ryanlayer/ssshtest

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run test_in_stdout python -c "print ('example assert_in_stdout success')"
assert_in_stdout "example"

run test_europa python examples/expand_europa.py
assert_no_stderr
assert_in_stdout "SparseEfficiencyWarning: Comparing sparse matrices using == is inefficient, try using != instead.
  x, y, x_list, y_list = netExp(R, P, x0, b)"
assert_exit_code 0
# assert_equal 'Data1.csv' $( ls 'Data1.csv' )
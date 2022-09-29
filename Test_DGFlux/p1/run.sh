#!/bin/bash
rm flow*dat
#../../src/dgsand input.dgsand1 >& log
../../src/dgsand input.dgsand1  input.dgsand2 >& log

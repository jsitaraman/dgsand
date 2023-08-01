#!/bin/bash

valgrind --leak-check=full \
         --show-leak-kinds=all \
         --track-origins=yes \
         --log-file=valgrind-out.txt \
         ../src/dgsand input.dgsand1  input.dgsand2 0.0 1.5
#         --verbose \

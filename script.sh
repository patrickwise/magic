#!/bin/bash
# This script is responsible for generating the parameter files for the basis sets of intrest
# and then handling the appropriate course of action.


function write_basis_set_file {
#opt: name maxn D b reserved_vars n_pol n_exp scaling_code data_file params_file der
echo "$2" > $1 || echo "Could not create file $1" && exit #maxn
echo "$3 $4" >> $1 #D b
echo "$5 $6 $7 $8" >> $1 #vars n_pol n_exp scaling
echo "$9 $10" >> $1 #data and param files
echo "$11 $12" >> $1 # deriv accuracy and deriv step size
echo "$13 $14" >> $1 # box size (xmin, xmax)
}

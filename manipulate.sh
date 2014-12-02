#!/bin/bash
#Manipulate a datafile to include the energy state and scaling and call plotter to plot it.
#Argument is the logfile, datafile.

#TODO: die if no argument.

energies=$(sed -n 6p $1 | sed 's/\t/ /g')



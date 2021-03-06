SYNOPSIS
=======================================


INSTALLATION
=======================================
Run:
``` bash
    make setup
```
in order to pull in dependency libraries.

Run:
``` bash
    make all
```
in order to build binary.


MOTIVATION
=======================================


RUN
=======================================
There are three modes to run in: c - compute a basis set indicated by
the run.file, r - run a ruggedness test with hard coded parameters, z -
run an experimental calculation that is hard coded into the binary.

The format of the run.file is:

|C Primitive      | Description                                                         |
|:--------------- |:------------------------------------------------------------------- |
|%s               | name of basis set                                                   |
|%zu              | number of excitation states                                         |
|%s               | name of nlopt algorithm to use for minimization                     |
|%lf %lf %lf      | D, b, mu; D and b of the Morse potential and mu is the reduced mass |
|%zu %lf          | accuracy and the step size of derivatives                           |
|%lf %lf %lf %lf  | integration parameters: xmin, xmax, epsrel and zero tolerance (s)   |
|%zu %zu %zu %c   | parameters of basis set function: growth, npol, nexp and scaling    |

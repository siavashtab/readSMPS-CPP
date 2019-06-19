#  readSMPS

[![Source](https://github.com/siavashtab/2slpEval/tree/master/2slpEval)](https://github.com/siavashtab/2slpEval/tree/master/2slpEval)
[![license](https://github.com/siavashtab/2slpEval/tree/master/LICENSE)](https://github.com/siavashtab/2slpEval/tree/master/LICENSE)

reading SMPS format files for two-stage stochastic programs
------------------

##  Author: Siavash Tabrizian -- stabrizian@gmail.com 

------------------

## Description:

readSMPS: A package for reading and saving the information of two stage stochastic programs from SMPS files

-  this package store the problem information in CPLEX Concert Technology format
-  the distribution of random variables are discrete (even if it is not discrete it can be turned to a discrete distribution)
-  randomness is on the right handside of the recourse function
-  COR file contains the mean value problem which its optimal value yields a lower bound

------------------

## Dependencies:

- This package needs the CPLEX(Concert Technology) solver for optimization parts. 

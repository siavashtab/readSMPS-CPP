#  readSMPS

[![Source](https://github.com/siavashtab/2slpEval/tree/master/2slpEval)](https://github.com/siavashtab/2slpEval/tree/master/2slpEval)
[![license](https://github.com/siavashtab/2slpEval/tree/master/LICENSE)](https://github.com/siavashtab/2slpEval/tree/master/LICENSE)

reading SMPS format files for two-stage stochastic programs

------------------

##  Author: Siavash Tabrizian -- stabrizian@gmail.com 

------------------

## Intoduction:

readSMPS: A package for reading and saving the information of two stage stochastic programs from SMPS files

It seems that there are implementations for reading SMPS files which are either written in other 
programming languages (C or Julia) or they do not provide a suitable data structures for L-shaped 
based algorithm (based on CPLEX solver). Moreover, it would be crucial to facilitate the problem that can handle sampling 
techniques for stochastic programming which can be done in this code.

-  this package store the problem information in CPLEX Concert Technology format
-  the distribution of random variables are discrete (even if it is not discrete it can be turned to a discrete distribution)
-  randomness is on the right handside of the recourse function
-  COR file contains the mean value problem which its optimal value yields a lower bound

------------------

## Dependencies:

- This package needs the CPLEX (Concert Technology) solver for optimization parts, 
  and the problem will be created based on the  CPLEX (Concert Technology) objects.

------------------

## Output

- this program will read the SMPS files, decompose the problem
  into a master and subproblem based on the CPLEX (Concert Technology) objects. 
  The data structure in which problems will be stored are all defined in prob_struct.h

------------------

## SMPS 

- SMPS is an efficient and effective way of describing stochastic programs. 

  It contains the follwoing files:
  
  1 - _.cor:
  
     This file is the core of the formulation which is basically derived from the 
	 mean value problem. It's format is very similar to _.mps files.
	 
  2 - _.tim
    
	This file contains the information of stages. The name of rows and columns 
	associated with each stage
	
  3 - _.sto
  
    This file contains the information about the random variables, and their distribution.
	Also, it shows the place in which the random variables appear in the second stage problem

-------------------

## Description

The program starts with asking the instance name and the number of samples
~~~~
>> Inter Test Instance Name:
>> Sample Size:
~~~~
After that it will produce the decomposed problems which is suitable for 
L-shaped type of algorithms


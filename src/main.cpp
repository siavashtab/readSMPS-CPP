/*****************************************************************************************\
**
** This code is written for reading SMPS format files and saving it in terms of the 
** CPLEX(Concert Technology) objects ( structures are defined in prob_struct.h)
**
**
**
**
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**           https://github.com/siavashtab
**   Copyright (c) 2019. All rights reserved.
\******************************************************************************************/

#include"ProbPrep.h"

void 
main()
{

	std::cout << "Inter Test Instance Name:";
	std::string instance;
	std::cin >> instance;

	std::cout << "Sample Size:";
	int initial_scen;
	std::cin >> initial_scen;

	//Initial the problem instance
	ProbPrep problem;
	problem.initialize(instance, initial_scen);


}

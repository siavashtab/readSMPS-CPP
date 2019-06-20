/*****************************************************************************************\
**
** This code is written for reading SMPS format files and solving them using ACS sampling
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**           https://github.com/siavashtab
**
\******************************************************************************************/

#include"ProbPrep.h"

void 
main()
{

	std::cout << "Inter Test Instance Name:";
	std::string instance;
	std::cin >> instance;

	//Initial the problem instance
	ProbPrep* problem;
	problem = new ProbPrep();
	problem->initialize(instance, initial_scen);

	problem->~ProbPrep();
	delete problem;


}
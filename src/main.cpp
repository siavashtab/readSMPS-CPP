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

#include"SAA.h"

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

	//Initial the SAA procedure
	SAA* saa;
	saa = new SAA();
	std::vector<SAA_output_struct> output_vec;
	saa->Initialize(problem);
	output_vec.push_back(saa->AlgOutput);
	saa->reset_SAA();

	problem->~ProbPrep();
	delete problem;

	if (Experiments == 1)
	{
		output_vec.clear();
		for (int expr = 0; expr < Expr_num; expr++)
		{
			ProbPrep* problem;
			problem = new ProbPrep();

			if (true) printf("expr: %d \n", expr);

			problem->initialize(instance, initial_scen + expr * initial_scen);
			saa->Initialize(problem);
			output_vec.push_back(saa->AlgOutput);
			saa->reset_SAA();

			problem->~ProbPrep();
			delete problem;

		}
	}

	delete saa;

	std::ofstream myfile(".\\output\\" + instance + "results.txt");

	myfile << "\n\n";
	myfile << "Ini scen,LB             ,UB            ,scen          ,pess Gap       ,opt time  \n";
	myfile << "----------------------------------------------------------------------------------------------\n\n";
	for (int i = 0; i < output_vec.size(); i++)
	{
		char buffer[50];
		sprintf(buffer, " %d , %0.2f +- %0.2f , %0.2f +- %0.2f , %0.2f +- %0.2f , %0.3f , %0.3f", initial_scen + i*initial_scen, 
			    output_vec[i].LB_mean, output_vec[i].LB_CI_std, output_vec[i].UB_mean, output_vec[i].UB_CI_std,
			output_vec[i].mean_scen_num, output_vec[i].std_scen_num, output_vec[i].Pessim_gap,
			output_vec[i].opt_time);
		myfile << buffer ;
		myfile << "\n\n";
	}
	myfile << "----------------------------------------------------------------------------------------------";
	
	system("pause");

}
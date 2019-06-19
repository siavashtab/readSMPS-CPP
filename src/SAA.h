
/*****************************************************************************************\
**
** SAA Class
**
** This file contains the routines needed to solve SP problems using SAA
**
** void optimize(ProbPrep&  problem);                             for optimization in 
**                                                                           each replication
** int eval(ProbPrep& problem, Sampling& sampl_ing, int siter);   for evaluation of the answer
**                                                                    obtained in each replication
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
\******************************************************************************************/

#ifndef SAA_h
#define SAA_h

#include "LShaped_Algo.h"


class SAA {  //Class for replicating the SP algorithms

public:

	SAA();
	void reset_SAA();

	int N;
	int saa_rep;
	int part_num;
	int ub_iter;
	int exp_num;

	double stdes;
	void optimize(ProbPrep*  problem);        // Main SAA Loop
	void Initialize(ProbPrep*  problem);
	void write_legend(std::ofstream& file, std::string filename);        // Print Solution
	std::vector<SOL_str>  x;

	double rep_LB_std;
	double rep_UB_std;
	double rep_LB_mean;
	double rep_UB_mean;
	double comp_UB_mean;
	double comp_UB_std;

	SAA_output_struct AlgOutput;

private:

	int eval(ProbPrep* problem, Sampling* sampl_ing, int siter);                       //Evaluation function 
	double opt_time;

	std::vector<double> rep_LB;
	std::vector<int> rep_partitions;
	std::vector<int> rep_Scen;
	std::vector<double> xhat;
	std::vector<double> rep_UB;

	std::vector<Partitions> parts;
	vec_d   UB_iter_firstIter;   //Algo upperbound of the first stage across replications
};



#endif //!SAA_h
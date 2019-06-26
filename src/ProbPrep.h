/*****************************************************************************************\
**
** ProbPrep Class (ProbPrep.h)
**
** This file contains the routines needed to read SMPS format files
**
**
**
**
**
** History:
**   March 2019 - Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
\******************************************************************************************/

#ifndef ProbPrep_h
#define ProbPrep_h


#include"Solver_CPLEX.h"
#include"RV_structs.h"
#include"commons.h"
//#include"CAP.h"

class ProbPrep
{

public:
	
	ProbPrep();
	~ProbPrep();

	SPProb_INFO * SPprobINFO;

	void initialize(std::string filename, int scen_num);
	void Probs_from_SMPS();     //create Probs by reading SMPS files

	Solver_CPLEX * solver;    //vector of solvers - everytime we need a solver we can create one

	//Mean Prob INFOS
	Prob * mean_prob;
	//***********************************************************

	//Created Master Probs
	Prob * master_prob;
	std::vector<SOL_str> x_reg;
	void         add_master_obj(Prob& prob, Solver_CPLEX& solver);
	void         add_surrogate_master_vars(Prob& prob, Solver_CPLEX& solver);
	//***********************************************************

	//Prob ph_sub_prob; 

	//Created Sub Probs
	std::vector<Prob> * stage_sub_prob;
	//***********************************************************


private:

	void printHeader();

	//Functions for reading the SMPS file format
	void read_COR();  //read the core file
	void read_TIM();  //read the time file
	void read_STOC(); //read the stoc file
	void read_STOC_INDEP(); //read the stoc file with dindep discrete
	//***********************************************************



	//Functions of generating master problems and its subroutins
	Prob         create_master_prob();              //create master prob
	void         add_master_vars(Prob& prob, Solver_CPLEX& solver);
	void         add_master_rngs(Prob& prob, Solver_CPLEX& solver);
	void         print_master_prob(Prob& prob);
	//***********************************************************


	//Functions of generating sub problems and its subroutins
	Prob         create_sub_probs(int stage);       //create stages sub probs
	void         add_sub_vars(Prob& prob, Solver_CPLEX& solver, int start, int end);
	void         add_sub_rngs(Prob& prob, Solver_CPLEX& solver);
	void         add_sub_obj(Prob& prob, Solver_CPLEX& solver);
	void         print_sub_prob(Prob& prob, int stage);
	void         fill_out_rhs(Prob& prob, Solver_CPLEX& solver);
	//***********************************************************

	void newProb(Prob &tmpProb);
	void newSPparam();
	void newRV(RV_info &rv);

};



#endif // !ProbPrep_h

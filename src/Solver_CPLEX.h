/*****************************************************************************************\
**
** Solver_CPLEX Class (Solver_CPLEX.h)
**
** This file contains the routines needed to use the optimization solver(CPLEX)
**
**
**
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
\******************************************************************************************/

#ifndef Solver_CPLEX_h
#define Solver_CPLEX_h

#include"prob_struct.h"
#include"Config.h"

class Solver_CPLEX {

public:
	
	std::string name;

	void Create_Var(Var_struct& var);	
	void Create_RNG(RNG_struct& rng);	
	void Create_OBJ(OBJ_struct& obj);	
	void Create_Prob(Prob& prob);
	void Clear_Prob(Prob& prob, std::string type);

	void open_prob(Prob& prob);
	void open_solver();
	void open_rng(IloRangeArray& rng, int iter);
	void end_solver();
	
	std::string getProbtype(Prob& prob);
	void setDefault(Prob& prob);
	

	IloExpr set_QP_obj(Prob& probl, std::vector<SOL_str> xhat);
	void get_Linear_obj_coeffs(Prob& prob);
	void get_Linear_rng_coeffs(Prob& prob);

	//Optimization related Functions Duals, Rays, ...
	IloNum GetDual(IloRange& rng, IloCplex& cplex);
	IloNum GetDual_Var_Bound(IloNum val, IloNumVar var);
	//*************************************************

	//---Usual Functions: Name, ...
	void ImportModel(Prob& prob);
	const char *getName(IloExtractable &any_object);
	int getSize(IloArray<IloExtractable> &any_object);
	void add_to_array(IloArray<IloExtractable> &any_object_array, IloExtractable &any_object);
	bool Solve_Prob(Prob& prob, bool isEval);
	void Print_Prob(Prob& prob); //(.lp)LP - (.mps)MPS
	void set_rhs_var(Prob& prob, std::vector<SOL_str>& var, std::vector<SOL_str>& rhs, std::vector<IloNum>& Cx);
	void set_rhs(Prob& prob, std::vector<SOL_str> rhs);
	int  rng_type(IloRange& rng);
	void decompose_range(IloRangeArray& mean_rng, Prob& prob, IloNumVarArray vars, int start, int end);
	void set_rhs_val(IloRange& rng, IloNum val);
	IloNum get_rhs(IloRange& rng);
	vec2_d getRHSSA_rand(Prob& prob);
	vec2_d getRHSSA_rand(Prob& prob, vec_d& lowe, vec_d& upper);
	void multip_d_toRange(double mult, IloRange& rnge);

	//creating cplex objects
	IloNumVar Create_Var_explicit(IloNum lb, IloNum ub, int type, const char *name);
	IloRange Create_RNG_explicit(IloNum lb, IloNum ub, int type, IloExpr expr, const char *name);
	IloObjective Create_OBJ_explicit(IloObjective::Sense sense, IloExpr expr, const char *name);
	//*****************************
	//*****************************

	//creating cutting plane related funcitons
	void AddCutToModel(Prob& prob);
	//***************************

	//******************************
	void Error_Report(Prob& prob, std::string Error_report);
	//*******************************

private:

	IloEnv env;


};

#endif // !Solver_h

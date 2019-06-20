

#include"Header.h"

//Sparse Vector Structure
struct sparse {
	int type;         //(0)int - (1) float - (2) double
	int pos_size; //index size
	vec_i pos;    //position of the positive values
	double dvalues; //values correspond to those positions
	int    ivalues; //values correspond to those positions
	float  fvalues; //values correspond to those positions
	const char *name;
	std::string id;
};
//**************************************************************************

struct Coeff_Sparse 
{
	int row;
	int col;
	const char *row_name;
	const char *col_name;
	IloNum val;

	bool operator<(const Coeff_Sparse& a) const
	{
		return col < a.col;
	}

};
//**************************************************************************

typedef std::vector<sparse> vec_sprs;

//variable structure
struct Var_struct {
	int type;        //0(bool) - 1(num) - 2(int)
	int pos_size;    //index size (for example pos_size = 3 like x_0_1_2 )
	vec_i pos;       //indexes
	IloNumVar x;
	const char *name;
	int lb;
	int ub;
	std::string id;
	IloNum obj_coef;
	std::vector<IloNum> rng_coefs;
	std::vector<std::string> rng_names;
};
//**************************************************************************

struct prev_Stage {
	vec_sprs coef;
	vec_sprs var_value;
	const char *var_name;
	std::string id;
	int stage;
};
//**************************************************************************

typedef std::vector<prev_Stage> PREV_STAGES; //if there are multiple variable values

//Solution structure when an optimization solver will be solved
struct SOL_str {
	int row_num;
	int col_num;
	const char *row;
	const char *col;
	IloNum value;
	IloNumVar *var;
	bool isRandom;
	int randomVar_indx;
	vec_d key;
};
//**************************************************************************

//constraints structure
struct RNG_struct {
	int type;        //0(<=) - 1(>=) - 3(==)
	int pos_size;    //index size (for example pos_size = 3 like x_0_1_2 )
	vec_i pos;       //indexes
	IloNum lb;
	IloNum ub;
	IloExpr expr;
	vec_sprs coef;
	PREV_STAGES prev;  //inf of the previous stage var
	IloRange rng;
	const char *name;
	std::string id;
	bool israndom;
	int rand_type;    //(0)rhs - (1) coeffs
};
//**************************************************************************

//objective function structure
struct OBJ_struct {
	int type;  //0(Max) - Min(1)
	IloNumArray coef;
	IloExpr     expr;
	IloNum      constant;
	IloObjective obj;
	const char *name;
	std::string id;
};
//**************************************************************************

//**************************************************************************
typedef std::vector<Var_struct> vec_var;
typedef std::vector<vec_var> vec_var2;
typedef std::vector<vec_var2> vec_var3;
typedef std::vector<vec_var3> vec_var4;
typedef std::vector<vec_var4> vec_var5;
typedef std::vector<RNG_struct> vec_rng;
//**************************************************************************

//**************************************************************************
struct var_vec {
	int dim;
	std::string id;
	vec_var var1;
	vec_var2 var2;
	vec_var3 var3;
	vec_var4 var4;
	vec_var5 var5;
};
//**************************************************************************

//Problem structure
struct Prob {
	IloEnv   * env;
	IloModel * model;
	IloCplex * cplex;
	bool isRaw;        //the raw prob is the one without seperated structures 
	                   //like rngs and vars and it just stored as a whole
	const char *name;
	std::string * strType;
	//--------------------------------
	IloRangeArray * range_raw;
	IloNumVarArray * vars_raw;
	bool has_surrogate;
	IloNumVarArray * surro_vars_raw;
	vec2_i       * surro_idx;           //index of the surrogate variables for adaptive sampling
	bool has_R;
	IloNumVarArray * R_vars_raw;
	IloObjective    * obj_raw;
	//---------------------------------

	//QP obj parameters----------------
	double lambda;
	double sigma;
	bool   update_xstar;
	//---------------------------------

	//Other Cuts ----------------------
	IloRangeArray * LShaped_opt;
	IloRangeArray * LShaped_opt_iter;
	IloRangeArray * LShaped_feas;
	IloRangeArray * master_cuts;
	std::vector<IloRange> * other_cuts;
	bool opt;
	bool feas;
	//---------------------------------

	const int * num_var;
	const int * num_rng;

	std::vector<Coeff_Sparse> * obj_coef_raw;
	std::vector<std::vector<Coeff_Sparse>> * rng_coefs_raw;

	//Previous Stages INFOs
	std::vector<std::vector<Coeff_Sparse>> * prev_rng_coefs_raw;
	std::vector<IloNumVarArray> * prev_vars_raw;

	int type; //0(LP) - IP(1) - MIP(2)
	std::string * probType;

	bool random;      //is it a random problem
	int type_random;  //(0)rhs - (1)obj - (3)coeffs - (4)rhsobj - ...
	std::string * id;

	//solution INFO------------------
	double zstar;
	double zstar_wout_surro;
	std::vector<SOL_str>  * sol;
	std::vector<SOL_str>  * xhat;  //if the answer in regularized around xhat
	bool             isReg; //is a regularized problem
	std::vector<SOL_str>  * sol_surrogate;
	std::vector<SOL_str>  * sol_R;
	std::vector<SOL_str>  * duals;
	vec_d            * dual_tot;// duals
	vec_d            * dual_R;  // dual of the random parts for random problems
	vec_d            * rho;     // r(w) - Cx
	vec_d            * Cx;      // Cx
	vec2_d           * beta;    // beta = C
	vec_d            * beta_sum;   // beta = sum(r in ranges, C(r))
	vec_d            * r_w;        // r(w)
	std::vector<SOL_str>  * opt_cut_duals;
	std::vector<SOL_str>  * feas_cut_duals;
	std::vector<SOL_str>  * rhs;
	//---------------------------------

	//print name
	std::string model_name;
};
//**************************************************************************

//Random Variable Structure
struct RV_info
{
	int rv_type;    //(0)RHS - (1)Cost - (2)Coeff
	int rv_row_num; //row number 
	int rv_col_num; //colum number 
	const std::string * ColName;
	const std::string * RowName;

	const char * name;
	std::string * id;

	//probability distribution info
	vec_d * val;
	vec_d * prob;

	//normal distribution
	double mean;
	double std;

	//uniform dist
	double l;
	double u;

	//CDF info
	std::map<double, double> * CDF;  //map<val, cumulative prob>

};
//**************************************************************************

//General information of the stohastic programs 
struct SPProb_INFO
{
	const int * num_var;
	const int * num_rngs;
	const int * rv_num;
	const int * stage_num;

	const int * SAA_iter;
	const int * initial_scen_num;
	const int * dynamic_scen_num;

	std::string file_name;
	std::string test_instance;
	std::string dirr_input;
	std::string dirr_algo;
	std::string dirr_model;
	std::string dirr_saa;
	std::string algorithm;
	std::string output_model;
	std::string output_algo;
	std::string output_saa;
	std::string model_type;
	std::string algo_type;

	//Time info:
	std::vector<std::vector<std::string>> * TIME_info;
	std::vector<int>            * TIME_col_idx;
	std::vector<int>            * TIME_row_idx;

	//StOC INFO
	std::vector<RV_info>  * RVs;
	std::string * STOC_TYPE;

};
//**************************************************************************

/*****************************************************************************************\
**
** ACS based Lshaped  Class
**
** This file contains the routines needed to solve a two-stage stochastic program
** using ACS
** 
**
**
**
**
** void optimize(ProbPrep& problem, Sampling& sampl_ing, vector<Partitions>& parts, int siter);
** void solve_sub(ProbPrep& problem, Sampling& sampl_ing, vector<Partitions>& parts);
** void solve_master(ProbPrep& problem);
**
**
** History:
**   Author: March 2019 - Siavash Tabrizian   stabrizian@gmail.com - stabrizian@smu.edu
**
**
\******************************************************************************************/

#include "Point_Estim.h"

struct LShaped_Cut {
	const char* Name;
	
	std::string cut_type;
	bool   isNew;            //if the cut is based on the new dual

	int num_scen;            //number of scen constructed the cut
	int  scen;             //vector of scenarios formed this cut

	double UB;      //upperbound associated with the cut
	double LB;      //lower bound associated with the cut


	double  dual_cut;  //dual of the cut in the master problem
	bool    isActive;  //is the cut active in optimality

	int iter;          //current iteration

	double Alpha_s;               //Intercept Info of the current iteration
	std::vector<double> Beta_s;        //Subgradient Info of the current iteration

	IloRange  cut;
	int cut_num;

	double weight;
	double multp;
};

struct LShaped_Info
{
	std::string type;

	std::vector<SOL_str> xFirst;
	std::vector<SOL_str> x_star;

	vec2_d opt_pi;                    //the optimal pis for each scenario
	std::vector<bool>   cut_status;        //if the cut is in the master
	vec_i          scens_to_remove;

	float sigma;
	double LB;
	double UB;
	double GAP;
	bool Algo_isStable;
	int max_iter;
	int min_iter;

	vec_d UB_now;        //vector of the upper bounds without weight
	vec_d UB_w_iter;     //vector of the upper bounds with weight
	double  UB_std_iter;
	double  UB_iter;
	vec_d LB_iter;       //vector of the lower bounds with weight
	std::vector<double>  GAP_it;
	vec_i  elite_samples;  //samples have found multiple partitions

	//adaptive Bendres Infos
	//Within Variability stability parameters
	double SSW;                      //sum of square within cluster
	vec_d   SSW_vec;                 //sum of square within cluster
	double  SSW_var;                 //variablity of sum of square within cluster

};

struct Cluster_Sample {
	vec2_i scens_iter;
	vec_i scens;
	vec_d  pi;
	int    size;
	double weight;
};


class LShaped_Algo {

public:

	LShaped_Algo();
	~LShaped_Algo();

	LShaped_Info LShapedINFO;

	double Alpha;
	std::vector<double> Beta;

	std::vector<double> Alpha_s;
	std::vector<std::vector<double>> Beta_s;

	void optimize(ProbPrep* problem, Sampling* sampl_ing, std::vector<Partitions>& parts, int siter);

	std::vector<std::vector<IloNum>> pi_collect;      //All duals 

	int  cut_num;
	double opt_time; //optimization time
	int cumul_tot_scen;

private:

	void Initialize(ProbPrep* problem);

	int iter;
	int nScenario;                           //number of initial scenarios

	std::vector<std::string> cut_type_it;
	std::vector<IloNum> pi;                      // current dual vector
	

	void solve_sub(ProbPrep* problem, Sampling* sampl_ing, std::vector<Partitions>& parts);
	void solve_master(ProbPrep* problem);
	void Point_Estimation_Analysis(ProbPrep* problem, Sampling* sampl_ing);

	//LShaped cut preparation
	void put_opt_in_range(ProbPrep* problem);
	void create_cuts(ProbPrep* problem);
	//-------------------------------------------------------------

	int num_partition_iter;              //number of clusters at the end of the iteration

	std::vector<std::vector<LShaped_Cut>> LShaped_cut_tot;     //LShaped cluster cuts for all iterations
	std::vector<LShaped_Cut> LShaped_cut;                 //LShaped cluster cuts for the current iteration

	bool stopping_check();                           //checking the stopping criteria for Adaptive SD

	//Type of generated Cuts:
	bool Optimality;
	bool feasibility;
	bool new_samples_generated;

	//Cluster Sampling Funcitons
	std::vector<Cluster_Sample>  samp_clust;
	int    tot_scen; //total number of scenarios
	int    tot_clust;
	double tot_norm;
	vec_d  scen_pi_norm;
	vec_i  scen_select;

	//for cluster based ACS we need to kep track of the cluster sizes
	int size_prev;
	int size_next;

	//Point Estimation class
	Point_Estim* p_estim;
	
};
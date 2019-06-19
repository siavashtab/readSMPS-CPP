/*****************************************************************************************\
**
** Point_Estimate Class
**
**
**
**
** This file contains the routines needed to estimate the recourse function
** for a given first stage variables (x)
**
**
**
** void monte(ProbPrep& problem, Sampling& sampl_ing); for monte carlo sampling method
** void lhs(ProbPrep& problem, Sampling& sampl_ing);   for LHS method
** void acs(ProbPrep& problem, Sampling& sampl_ing);   for ACS method
**
**
**
**
** History:
**   Author: March 2019 - Siavash Tabrizian   stabrizian@gmail.com - stabrizian@smu.edu
**
\******************************************************************************************/



#ifndef Point_Estim_h
#define Point_Estim_h

#pragma once

#include "Sampling.h"
#include "utils.h"

//Scenario structure *************************************
struct Scen_struct {

	double est;
	vec_d  *samp;
	bool   active;
	double alpha;
	vec_d  *beta;
	vec_d  *z1;

};
//*******************************************************

//Scenario Cluster Structure ****************************
struct Scen_cluster {

	double estim;
	double estim_R;   //estimate the random part  = r(w)^T * pi_R 
	double estim_D;   //estimate the det part  =    estim - estim_R
	double estim_lower;
	double estim_upper;
	vec_d  estim_extrm;
	std::vector<Scen_struct> *scens; //collection of scenarios in each cluster
	int   size; //number of scenarios in each cluster
	vec_i scen_idx; 
	bool  active;
	double sum_freq;

	//samples in cluster
	vec2_d *samples;
	vec2_d *samples_extrm;
	vec_d  *lower;
	vec_d  *upper;
	vec_i  *size_of_samps;
	vec2_d *z1;     //zero one interval of each random variable
	double SSW;    //sum squared error within cluster
	vec_d  SSW_it;    //sum squared error within cluster


	//Attribute of the cluster:
	vec_d  *pi;           //pi of the random part
	int maxr_;            //rv with the highest pi value
	vec_d  *pi_tot;       //pi 
	vec2_d *SA_vals;      //sensitivity range
	vec2_d *SA_z1;        //zero one interval of each random variable of SA
	int    max_pi;
	bool   new_pi;
	int pos;              //position of the dual in the pool
	int clust_size;
	vec_d r_Cx;
	vec_d Cx;
	int pi_num;

	//cut info
	double  alpha_estim;    // estimate of alpha = pi^T * r(w)
	vec_d   *beta_estim;    // estimate of beta = pi^T * C
	double  weight;         // 1/(cluster size * number of clusters)

};
//*******************************************************


class  Point_Estim {

public:

	Point_Estim();
	~Point_Estim();

	void initialize(int sampnum, int rep);
	void end();

	std::vector<SOL_str>  *x;


	//Output from differnt sampling techniques
	void monte(ProbPrep* problem, Sampling* sampl_ing);
	void lhs(ProbPrep* problem, Sampling* sampl_ing);
	void acs(ProbPrep* problem, Sampling* sampl_ing);
	//*******************************************************

	//explicit data structures
	double  estimate;
	vec_d   *estimate_c;    //cluster estimate
	vec2_i  *scen_idx;      //scenario numbers of each cluster
	double  estimate_std; 
	int     scen_num; // scen number which sequentially increases
	int     n_c_tot;  //tot number of scenarios
	double  mean_MSW; //mean of MSW of each cluster
	double  R_a;
	double  tot_size; //total cluster size
	int dual_num;
	vec2_d  *samples;
	//*******************************************************

	std::vector<Scen_cluster>  *cluster_samples;

private:

	//Parameters
	int     samp_num;
	vec2_d  *dual_pool;
	vec_i   *dual_pool_num;
	vec_i   *dual_pool_clust;
	int     initial_samp;
	vec_i   *samples_clust; //which cluster samples belong
	vec_d   *scen_estim;
	bool    NoClust; //if no cluster can be discovered in seq function
	vec_d   clustNumIt; //number of clusters in each iter of seq
	vec_d   MSW_it;
	bool    isACS; //sampling is ACS

	//Seed
	long long seed;

	//functions needed for ACS
	void create_neighbor(ProbPrep* problem, Sampling* sampl_ing, int clust, int num);
	void create_neighbor_None(ProbPrep* problem, Sampling* sampl_ing, int clust, int num);
	void create_neighbor_all(ProbPrep* problem, int clust);
	void create_neighbor_None_rmax(ProbPrep* problem, Sampling* sampl_ing, int clust, int num);
	void seq(ProbPrep* problem, Sampling* sampl_ing);
	void update_estimate();
	void update_estimate_seq();
	std::vector<Scen_struct> lhs2(ProbPrep* problem, Sampling* sampl_ing, vec2_d z1, int size);
	std::vector<Scen_struct> monte2(ProbPrep* problem, Sampling* sampl_ing, vec2_d z1, int size);
	bool is_in_SArange(vec2_d& SA, vec_d& z1);
	double getSize(Scen_cluster& clust);
	void   Req_neighb_samp(vec2_d& tot_samp, vec_d& orig_samp, int start, ProbPrep* problem, int clust);
	vec_d   neighb_samp(vec_d& orig_samp, ProbPrep* problem, int rv, int updwn);
	void    add_scen(ProbPrep* problem, vec_d& samp, int clust);
	void    add_scen2(ProbPrep* problem, vec_d& samp, int clust);
	void    create_z1(vec2_d& z1_main, vec_d& z1);
	void    update_z1(vec2_d& z1_main, vec_d& z1);
	//*******************************************************

	

};


#endif  //!Point_Estim
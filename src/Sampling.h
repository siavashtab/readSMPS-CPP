
/*****************************************************************************************\
**
** Sampling Class
**
** This file contains the routines needed to generate samples
** 
**
**
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
\******************************************************************************************/


#ifndef Sampling_h
#define Sampling_h

#include"ProbPrep.h"


class  Sampling
{

public:

	Sampling();
	~Sampling();


	std::vector<SSample> samples;                               //all samples structures
	std::vector<Clusters_SSample> clusterSamples;               //samples put in clusters
	std::vector<Partitions>       partitions;                   //collection of dual partitions

	void set_random_rhs(std::vector<SOL_str>&  rhs, int scen);
	void set_random_rhs_EVAL(std::vector<SOL_str>&  rhs, std::vector<RV_info>& rv, std::mt19937& randgen);
	
	//PE: Point Estimation - type:0(monte) 1(lhs) 2(adaptive)
	void set_random_rhs_PE(vec_d& emptsample, vec_d& z1, std::vector<SOL_str>&  rhs, std::vector<RV_info>& rv,
		                   std::mt19937& randgen, vec_d Ll, vec_d Ul);
	vec_d   real_in_range(RV_info& rv, vec_d& range);
	vec2_d  inRHSSA(std::vector<RV_info>& rv, vec2_d val);
	vec2_d  inRHSSA(std::vector<RV_info>& rv, vec2_d val, vec_d& lower, vec_d& upper);


	int initial_scen_num;

private:

	//generating a random realization from a given 0-1 interval
	IloNum RV_Gen(RV_info& rv, std::mt19937& randgen, vec_d& zero_one_l, bool adjust_interval);
	vec_d  rev_zero_one(RV_info& rv, vec_d val);
	vec_d  rev_zero_one(RV_info& rv, vec_d val, double& lower, double& upper);
	double  rev_zero_one_d(RV_info& rv, double val);
	IloNum RV_Gen_EVAL(RV_info& rv, std::mt19937& randgen);
	std::tuple<IloNum, double> RV_Gen_PE(RV_info& rv, std::mt19937& randgen, double Ll, double Ul);
	//**********************************************************

	double LHS_0_1_length;

	//sampling sub functions
	vec_d zero_one_LHS(double val);


	double sample_variablity(std::vector<SSample>& samples);
	
};


#endif // !Sampling_h
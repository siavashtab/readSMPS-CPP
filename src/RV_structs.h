
#include"Header.h"


struct SAMPLE 
{

	std::vector<IloNum> Samples;
	vec2_d LHS_inf;
	double Prob;

};

struct SSample
{
	int main_scen;

	bool   isMain;
	bool   isRef;
	int    cluster_sample_idx;  //index is scen struct
	int    ref_scen;    //reference scen
	vec_i  sub_scen;    //sub scens
	SAMPLE sample;      //SAMPLE info
	int    freq;    

	bool   isOnBorder; 
	bool   isMax;
	int    RV_num_border;

	int    new_parts_found;  //number of new partitions this scenario has found
	bool   isdeActive;  
	int    depth;       //How many heirarchies passed
	double weight;
	std::string Method;
	bool isAdaptSample;

	//Optimization Infos
	double UB;
};

struct Clusters_SSample
{
	int mainscen;                       //the main scenario
	int size;                           //number of scenarios in the cluster
	std::vector<SSample> samples;    //all sample info in the clusters
	double weight;

	double UB;                          //UB of the cluster
	vec_d  sub_UB;                      //all UB of scenarios in the cluster

	double within_cluster_square;
	int new_partition_num;              //number of new partions this sample found               
	bool NewPartitionsFound;
	bool Active;
	int SubScenDepth;

	double pi_var;
	double sample_var;

	double rpi;          //the term sum(r . pi - sum( r . pi ) )
	double Cpi;          //the term sum(C . pi . x - sum( C . pi . x ) )
};

struct Partitions
{
	double   pi_norm;
	vec_d    pi;
	vec_i    scens;           //scenarios inside this partition in the last iteration
	vec_i    max_extreme_scens;
	vec_i    min_extreme_scens;
	vec_i    iter;            //iterations that this partition found
	bool    isActive;         //is the pi is important in the LShaped cut
};
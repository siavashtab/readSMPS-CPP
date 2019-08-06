
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



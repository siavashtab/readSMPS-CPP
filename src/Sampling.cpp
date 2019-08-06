
/*****************************************************************************************\
**
** Sampling Class (Sampling.cpp)
**
** This file contains the subroutines needed to generate samples
**
**
**
**
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**   Copyright (c) 2019. All rights reserved.
\******************************************************************************************/

#include"Sampling.h"
#include"seeds.h"

Sampling::Sampling()
{
	std::cout << "--------------------------------------------------------------------------------" << std::endl;

	std::cout << "                               INITIAL SAMPLING                                 " << std::endl;

	std::cout << "--------------------------------------------------------------------------------" << std::endl;
}

Sampling::~Sampling()
{
	samples.clear();
	clusterSamples.clear();
	partitions.clear();
}


//Generating random number inside a zero one interval
IloNum Sampling::RV_Gen(RV_info& rv, std::mt19937& randgen, vec_d& zero_one_l, bool adjust_interval)
{
	std::uniform_real_distribution<double> dis(zero_one_l[0], zero_one_l[1]);

	double rnd = (double)dis(randgen);

	if (Display_Samples == 1) std::cout << rnd << " , ";

	std::map<double, double>::iterator it;

	for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
	{
		if (rnd <= it->second)
		{
			if(adjust_interval == true) zero_one_l = zero_one_LHS(it->first);
			return (IloNum)it->first;
			break;
		}
	}

	return (IloNum)it->first;
	
}

//generating samples inside the sensitivity analysis bounds
vec2_d  Sampling::inRHSSA(std::vector<RV_info>& rv, vec2_d val)
{
	vec2_d res;
	vec_d  res1;


	if (rv.size() != val.size())
	{
		printf("ERROR: SAmpling: max_inRHSSA: not equal sizes");
		system("pause");
	}
		

	for (int r = 0; r < rv.size(); r++)
	{
		vec_d zo_interval = rev_zero_one(rv[r], val[r]);
		//printf("rv :  %d , lower: %0.02f ,  upper: %0.02f  \n", r, zo_interval[0], zo_interval[1]);
		res1.push_back(std::max(0.0,zo_interval[0] - 0.0001));
		res1.push_back(std::min(1.0, zo_interval[1] + 0.0001));
		res.push_back(res1);
		res1.clear();
	}



	return res;
}

//generating samples inside the sensitivity analysis bounds
vec2_d  Sampling::inRHSSA(std::vector<RV_info>& rv, vec2_d val, vec_d& lower, vec_d& upper)
{
	vec2_d res;
	vec_d  res1;


	if (rv.size() != val.size())
	{
		printf("ERROR: Sampling: max_inRHSSA: not equal sizes \n");
		printf("val size: %d  -  rv size:  %d", val.size(), rv.size());
		system("pause");
	}


	for (int r = 0; r < rv.size(); r++)
	{
		vec_d zo_interval = rev_zero_one(rv[r], val[r], lower[r],upper[r]);
		if(false) printf("rv :  %d , lower: %0.02f ,  upper: %0.02f  \n", 
			                                     r, zo_interval[0], zo_interval[1]);
		res1.push_back(std::max(0.0, zo_interval[0]));
		res1.push_back(std::min(1.0, zo_interval[1]));
		res.push_back(res1);
		res1.clear();
	}



	return res;
}

//Get the zero one interval for a given values
vec_d Sampling::rev_zero_one(RV_info& rv, vec_d val)
{
	vec_d tmp;
	tmp.resize(2);
	tmp[0] = 0.0;
	tmp[1] = 1.0;

	std::map<double, double>::iterator it;

	for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
	{
		if (val[0] <= it->first)
		{
			tmp[0] = it->second;
			break;
		}
		else
		{
			tmp[0] = it->second;
		}
	}

	for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
	{
		if (val[1] >= it->first)
		{			
			tmp[1] = it->second;
		}
		else
		{
			tmp[1] = it->second;
			break;
		}
	}

	if (false)
	{
		printf("val0: %0.02f  ,  val1: %0.02f , z1_0: %0.06f , z1_1: %0.06f  \n", val[0], val[1], tmp[0], tmp[1]);
		int ss = getchar();
	}


	return tmp;

}

//Get the zero one interval for a given values
vec_d Sampling::rev_zero_one(RV_info& rv, vec_d val, double& lower, double& upper)
{
	vec_d tmp;
	tmp.resize(2);
	tmp[0] = 0.0;
	tmp[1] = 1.0;

	std::map<double, double>::iterator it;

	if (val[0] == val[1])
	{
		if (val[0] == 0)
		{
			tmp[0] = 0.0;
			tmp[1] = 0.0;
			lower = 0.0;
			upper = 0.0;
		}
		else
		{
			for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
			{
				if (val[0] >= it->first)
				{
					lower = it->first;
					upper = it->first;
					tmp[0] = it->second;
					tmp[1] = it->second;
				}
				else
				{
					break;
				}

			}
		}

	}
	else
	{
		for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
		{
			if (val[0] >= it->first)
			{
				lower = it->first;
				tmp[0] = it->second;			
			}
			else
			{
				break;
			}

		}

		for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
		{
			if (val[1] >= it->first)
			{
				upper = it->first;
				tmp[1] = it->second;
			}
			else
			{
				break;
			}
		}
	}



	return tmp;

}

//Get the zero one interval for a given values
double Sampling::rev_zero_one_d(RV_info& rv, double val)
{
	double tmp;
	tmp = 1.0;

	std::map<double, double>::iterator it;

	for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
	{
		if (val > it->first)
		{
			tmp = it->second;
			break;
		}
		else
		{
			tmp = it->second;
		}
	}

	return tmp;

}

IloNum Sampling::RV_Gen_EVAL(RV_info& rv, std::mt19937& randgen)
{
	std::uniform_real_distribution<double> dis(0, 1);
	double rnd = (double)dis(randgen);

	for (std::map<double, double>::iterator it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
	{
		//cout << it->second << ",";
		if (rnd <= it->second)
		{
			return (IloNum)it->first;
			break;
		}
	}

	printf("=");

	std::map<double, double>::iterator it = rv.CDF->end();
	return (IloNum)it->first;

}

std::tuple<IloNum, double> Sampling::RV_Gen_PE(RV_info& rv, std::mt19937& randgen, double Ll, double Ul)
{
	std::uniform_real_distribution<double> dis(Ll, Ul);
	double rnd = (double)dis(randgen);
	bool printing = false;

	if (printing == true) std::cout << Ll << " , ";

	std::map<double, double>::iterator it;
	
	for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
	{
		
		if (rnd <= it->second)
		{
			if (printing == true) std::cout << rnd << " , ";
			if (printing == true) std::cout << Ul << std::endl;
			break;
		}
	}
	
	if (printing == true) std::cout << it->first << " " << it->second << std::endl;

	return std::make_tuple((IloNum)it->first, it->second );

}


vec_d Sampling::zero_one_LHS(double val)
{
	double sum = 0;
	vec_d empt;
	int end = ceil(val / LHS_0_1_length);
	for (int i = 0; i < end + 1; i++)
	{
		sum += LHS_0_1_length;
		if (val <= sum)
		{
			empt.push_back(sum - LHS_0_1_length);
			empt.push_back(sum);
			return empt;
			break;
		}
	}
}



void Sampling::set_random_rhs(std::vector<SOL_str>&  rhs, int scen)
{
	int loc = 0;
	if(!(rhs.size() > 0)) std::cout << "ERROR: Set RHS: rhs is empty!" << std::endl;
	for (int i = 0; i < rhs.size(); i++)
	{
		//cout << i << " ";
		if (rhs[i].isRandom)
		{
			if (rhs[i].randomVar_indx != loc) std::cout << "ERROR: Set RHS: referencing the rhs struct is not right!" << std::endl;
			rhs[i].value = samples[scen].sample.Samples[loc];
			//cout << rhs[i].value << ",";
			loc++;
			if (loc > samples[scen].sample.Samples.size()) break;
		}
		//cout << endl;
	}
}

void Sampling::set_random_rhs_EVAL(std::vector<SOL_str>&  rhs, std::vector<RV_info>& rv, std::mt19937& randgen)
{
	for (int i = 0; i < rhs.size(); i++)
	{
		//cout << i;
		if (rhs[i].isRandom)
		{
			//cout << " Random, idx = " << rhs[i].randomVar_indx;
			rhs[i].value = RV_Gen_EVAL(rv[rhs[i].randomVar_indx], randgen);
			//cout << " val: " << rhs[i].value ;
		}
		//cout << endl;
	}
}


void Sampling::set_random_rhs_PE(vec_d& emptsample, vec_d& z1, std::vector<SOL_str>&  rhs, std::vector<RV_info>& rv,
	                             std::mt19937& randgen, vec_d Ll, vec_d Ul)
{
	int rv_idx = 0;
	
	bool printing = false;

	for (int i = 0; i < rhs.size(); i++)
	{

		if(printing == true) std::cout << i;

		if (rhs[i].isRandom == true)
		{
			if (printing == true)std::cout << " Random, idx = " << rhs[i].randomVar_indx;
			double z_one;
			if(false) std::cout << rhs[i].randomVar_indx << " " << Ll[rv_idx] << " " << Ul[rv_idx] << std::endl;
			std::tie(rhs[i].value, z_one) = RV_Gen_PE(rv[rhs[i].randomVar_indx], randgen, Ll[rv_idx], Ul[rv_idx]);
			z1.push_back(z_one);
			emptsample.push_back((double)rhs[i].value);
			rv_idx++;
			if (printing == true) std::cout << " val: " << rhs[i].value ;
		}
		if (printing == true) std::cout << std::endl;
	}
	
}


vec_d Sampling::real_in_range(RV_info& rv, vec_d& range)
{
	vec_d output;

	std::map<double, double>::iterator it;
	bool print = false;

	if (print == true)
	{
		printf("LB: %0.02f ,", range[0]);
	}

	//map<val, cumulative sum>
	for (it = rv.CDF->begin(); it != rv.CDF->end(); ++it)
	{

		if (range[0] <= it->first && range[1] >= it->first)
		{
			output.push_back(it->first);
			if (print == true)
			{
				printf(" , %0.02f ", it->first);
			}
		}
	}

	if (print == true)
	{
		printf(", UB:  %0.02f  \n", range[1]);
		int ss = getchar();
	}

	return output;
}

#ifndef HEADER_H
#define HEADER_H
#pragma warning(disable:4786)
#include <stdlib.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <vector>
#include <map>
#include <set>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <numeric>
#include <iomanip>
#include <ilcplex/ilocplex.h>
#define INFINITY 99999999999

typedef std::vector<IloBoolVar> vec_bool;
typedef std::vector<vec_bool> vec_2bool;

typedef std::vector<IloNumVar> vec_num;
typedef std::vector<vec_num> vec_2num;

typedef std::vector<IloIntVar> vec_int;
typedef std::vector<vec_int> vec_2int;

typedef std::vector<int> vec_i;
typedef std::vector<double> vec_d;
typedef std::vector<float> vec_f;
typedef std::vector<vec_i> vec2_i;
typedef std::vector<vec_d> vec2_d;
typedef std::vector<vec2_d> vec3_d;
typedef std::vector<vec_f> vec2_f;


typedef IloArray<IloNumVarArray>  IloNumVarArray2;
typedef IloArray<IloBoolVarArray>  IloBoolVarArray2;
typedef IloArray<IloNumVarArray2>  IloNumVarArray3;
typedef IloArray<IloRangeArray>  IloRangeArray2;
typedef IloArray<IloNumArray> TwoDMatrix;
typedef IloArray<IloExprArray> IloExprArray2;
typedef IloArray<IloExprArray2> IloExprArray3;
typedef IloArray<IloExprArray3> IloExprArray4;

struct SAA_output_struct
{
	double LB_mean;
	double LB_std;
	double LB_CI_std;

	double UB_mean;
	double UB_std;
	double UB_CI_std;

	double Pessim_gap;

	double opt_time;
	double tot_time;

	double mean_part_num;
	double std_part_num;

	double mean_scen_num;
	double std_scen_num;
};


#endif;


/*****************************************************************************************\
** utils.h
**
**
**
** Utility functions for working with std::vector
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**   Copyright (c) 2019. All rights reserved.
\******************************************************************************************/


#include "Header.h"
#include <functional>

namespace utils {
	double L2_norm(std::vector<IloNum> const &u);
	double L2_norm(std::vector<double> const &u);
	double avg_func(std::vector<double> const &vv);
	double variance_func(const std::vector<double>& vec);
	double variance_func(const std::vector<int>& vec);
	bool   inVector(const std::vector<double>& vec, const double& tmp);
	bool   inVector(const vec2_d& vec, const vec_d& tmp);
	bool   inVector(const vec2_d& vec, const vec_d& tmp, int& pos);
	bool   inVector(const std::vector<int>& vec, const int& tmp);
	std::vector<int>    maxVector_num(const std::vector<double>& vec, int num);
	std::vector<int>    minVector_num(const std::vector<double>& vec, int num);
	std::vector<int>    meanVector_num(const std::vector<double>& vec, int num);
	int    maxVector(const std::vector<double>& vec);
	int    minVector(const std::vector<double>& vec);
	int    maxVector(const std::vector<int>& vec);
	int    minVector(const std::vector<int>& vec);
	double    meanVector(const std::vector<double>& vec);
	double    meanVector(const std::vector<int>& vec);
	double    stdVector(const std::vector<double>& vec);
	double    stdVector(const std::vector<int>& vec);
	double    stdVector(const vec2_d& vec);
	int    findVector(const std::vector<double>& vec, const double& tmp);
	int    findVector(const vec2_d& vec, const vec_d& tmp);
	vec_d  element_w_prod(const vec_d& v, const vec_d& u);
	vec_d  const_prod(const vec_d& v, double w);
	int    findintVector(const std::vector<int>& vec, int& tmp);
	double    l2normOfminus_vec(const std::vector<double>& v, const std::vector<double>& u);
	double    square_vec_d(const std::vector<double>& v);
	double    square_vec_d(const vec2_d& v);
	bool falls_in_range(const vec2_d& range, const vec_d& val);
	void print_vector(const vec_i& v);
	void print_vector(const vec_d& v);
	void print_vector(const vec2_d& v);
	double sum_square(const vec_d& v);
	vec_d  get_Centre(const vec2_d& v);
	double L1_dist(const vec_d& v, const vec_d& u);
	vec_d  prod_components(const vec_d& v, const vec_d& u);
	int    max_diff(const vec_d& v, double& val);
	int    max_diff(const vec_d& v, const vec_d& u);
}


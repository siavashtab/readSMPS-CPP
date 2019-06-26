
/*****************************************************************************************\
** utils.h
**
**
**
** Utility functions for working with vectors
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
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
	std::vector<int>    maxVector_num(const std::vector<double> vec, int num);
	std::vector<int>    minVector_num(const std::vector<double> vec, int num);
	std::vector<int>    meanVector_num(const std::vector<double> vec, int num);
	int    maxVector(const std::vector<double>& vec);
	int    minVector(const std::vector<double>& vec);
	int    maxVector(const std::vector<int>& vec);
	int    minVector(const std::vector<int>& vec);
	int    findVector(const std::vector<double>& vec, const double& tmp);
	int    findVector(const vec2_d& vec, const vec_d& tmp);
	vec_d  element_w_prod(vec_d& v, vec_d& u);
	vec_d  const_prod(vec_d& v, double w);
	int    findintVector(std::vector<int>& vec, int& tmp);
	double    l2normOfminus_vec(std::vector<double>& v, std::vector<double>& u);
	double    mean_vec_d(std::vector<double>& v);
	double    std_vec_d(std::vector<double>& v);
	double    std_vec_d(vec2_d& v);
	double    std_vec_i(std::vector<int>& v);
	double    square_vec_d(std::vector<double>& v);
	double    square_vec_d(vec2_d& v);
	bool compare_f(vec_d pi1, vec_d pi2, vec_d right, vec2_d beta);
	bool falls_in_range(vec2_d& range, vec_d& val);
	void print_vectori(vec_i& v);
	void print_vectord(vec_d& v);
	void print_vectord(vec2_d& v);
	double sum_square(vec_d& v);
	vec_d  get_Centre(vec2_d& v);
	double L1_dist(vec_d& v, vec_d& u);
	vec_d  prod_components(vec_d& v, vec_d& u);
	int    max_diff(vec_d& v, double& val);
	int    max2_diff(vec_d& v, vec_d& u);
}


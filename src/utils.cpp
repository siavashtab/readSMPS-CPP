
/*****************************************************************************************\
** utils.cpp
**
**
**
** Utility functions for working with vectors
**
** History:
**   Author: Siavash Tabrizian   stabrizian@gmail.com stabrizian@smu.edu
**
\******************************************************************************************/


#include "utils.h"

/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////

double utils::L2_norm(const std::vector<double>& u)
{
	double accum = 0.;
	for (int i = 0; i < u.size(); ++i) {
		accum += u[i] * u[i];
	}
	double norm = sqrt(accum);
	return norm;
}


//Function for average
double utils::avg_func(const std::vector<double>& vv)
{
	double return_value = 0.0;
	int n = vv.size();

	for (int i = 0; i < n; i++)
	{
		return_value += vv[i];
	}

	return (return_value / n);
}
//****************End of average funtion****************

bool utils::falls_in_range(const vec2_d& range, const vec_d& val)
{
	bool pass = true;

	for (int i = 0; i < range.size(); i++)
	{
		if (val[i] >= range[i][1] || val[i] <= range[i][0])
		{
			pass = false;
			break;
		}
	}

	return pass;
}

//Function for variance
double utils::variance_func(const std::vector<double>& vec)
{

	double mean = avg_func(vec);
	int n = vec.size();
	double sum = 0.0;

	for (int j = 0; j < n; j++)
	{
		sum += pow((vec[j] - mean), 2);
	}

	double var = sum / ((double)(n - 1));

	return var;
}
//****************End of variance funtion****************

//Function for variance
double utils::variance_func(const std::vector<int>& vec)
{

	int summ = accumulate(vec.begin(), vec.end(),0);
	double mean = summ / (double)vec.size();
	int n = vec.size();
	double sum = 0.0;

	for (int j = 0; j < n; j++)
	{
		sum += pow((vec[j] - mean), 2);
	}

	double var = sum / ((double)(n - 1));

	return var;
}
//****************End of variance funtion****************

// check if it is in a vector function

bool utils::inVector(const std::vector<double>& vec, const double& tmp)
{

	bool isInV;

	isInV = (std::find(vec.begin(), vec.end(), tmp) != vec.end());

	return isInV;


}

//****************End of check if it is in a vector funtion****************

// check if it is in a vector 2d function

bool utils::inVector(const vec2_d& vec, const vec_d& tmp)
{

	bool isInV = false;

	for (vec2_d::const_iterator it = vec.begin(); it != vec.end(); it++) {
		if (*it == tmp) 
		{
		   isInV = true;
		}
	}

	return isInV;

}

//****************End of check if it is in a vector 2d funtion****************


//****************End of check if it is in a vector 2Num funtion****************

// check if it is in a 2vector function

bool utils::inVector(const vec2_d& vec, const vec_d& tmp, int& pos)
{


	bool isInV = false;

	for (vec2_d::const_iterator it = vec.begin(); it != vec.end(); it++) {
		if (*it == tmp)
		{
			isInV = true;
			pos = distance(vec.begin(), it);
			break;
		}
	}

	return isInV;

}

//****************End of check if it is in a 2vector funtion****************


// check if it is in a vector function

bool utils::inVector(const std::vector<int>& vec, const int& tmp)
{

	bool isInV;

	isInV = (std::find(vec.begin(), vec.end(), tmp) != vec.end());

	return isInV;


}

//****************End of check if it is in a vector funtion****************

// find the position of a tmp in a vector
int utils::findVector(const std::vector<double>& vec, const double& tmp)
{

	std::vector<double>::const_iterator locc;
	locc = std::find(vec.begin(), vec.end(), tmp);
	int loc = distance(vec.begin(), locc);

	return loc;

}
//****************End of find the position of a tmp in a vector funtion****************

int utils::findVector(const vec2_d& vec, const vec_d& tmp)
{
	vec2_d::const_iterator locc;
	locc = std::find(vec.begin(), vec.end(), tmp);
	int loc = distance(vec.begin(), locc);

	return loc;
}

// find the position of a tmp in a vector
int utils::findintVector(const std::vector<int>& vec, int& tmp)
{

	std::vector<int>::const_iterator locc;
	locc = std::find(vec.begin(), vec.end(), tmp);
	int loc = distance(vec.begin(), locc);

	return loc;

}

// find maximum value in a vector
int utils::maxVector(const std::vector<double>& vec)
{

	std::vector<double>::const_iterator result;
	result = std::max_element(vec.begin(), vec.end());
	int loc = std::distance(vec.begin(), result);

	return loc;
}
//****************End of find maximum value in a vectorfuntion****************

// find a number of maximum values in a vector
std::vector<int> utils::maxVector_num(const std::vector<double>& vec, int num)
{

	int N = vec.size();
	std::vector<double> tmp;
	tmp = vec;
	std::vector<int> result;
	sort(tmp.begin(), tmp.end());
	for (int i = N - 1; i > N - num - 1; i--){
		int loc = utils::findVector(vec, tmp[i]);
		result.push_back(loc);
	}

	return result;
}
//****************End of find maximum value in a vectorfuntion****************

// find a number of minimum values in a vector
std::vector<int> utils::minVector_num(const std::vector<double>& vec, int num)
{

	int N = vec.size();
	std::vector<double> tmp;
	tmp = vec;
	std::vector<int> result;
	sort(tmp.begin(), tmp.end());
	for (int l = 0; l < num; l++){
		int loc = utils::findVector(vec, tmp[l]);
		result.push_back(loc);
	}

	return result;
}
//****************End of find minimum value in a vectorfuntion****************
// find a number of minimum values in a vector
std::vector<int> utils::meanVector_num(const std::vector<double>& vec, int num)
{

	int N = vec.size();
	std::vector<double> tmp;
	tmp = vec;
	std::vector<int> result;
	sort(tmp.begin(), tmp.end());
	int ll = round(N / 2) - floor(num / 2);
	int uu = round(N / 2) + floor(num / 2) + 1;

	for (int l = ll; l < uu; l++){
		int loc = utils::findVector(vec, tmp[l]);
		result.push_back(loc);
	}

	return result;
}
//****************End of find minimum value in a vectorfuntion****************
// find minimum value in a vector
int utils::minVector(const std::vector<double>& vec)
{

	std::vector<double>::const_iterator result;
	result = std::min_element(vec.begin(), vec.end());
	int loc = std::distance(vec.begin(), result);

	return loc;
}
//****************End of find minimum value in a vectorfuntion****************

// find maximum value in a vector
int utils::maxVector(const std::vector<int>& vec)
{

	std::vector<int>::const_iterator result;
	result = std::max_element(vec.begin(), vec.end());
	int loc = std::distance(vec.begin(), result);

	return loc;
}
//****************End of find maximum value in a vectorfuntion****************

// find minimum value in a vector
int utils::minVector(const std::vector<int>& vec)
{

	std::vector<int>::const_iterator result;
	result = std::min_element(vec.begin(), vec.end());
	int loc = std::distance(vec.begin(), result);

	return loc;
}
//****************End of find minimum value in a vectorfuntion****************

// || v - u||^2
double utils::l2normOfminus_vec(const std::vector<double>& v, const std::vector<double>& u)
{

	double d_uv = 0;
	int v_size = v.size();

	for (int n = 0; n < v_size; n++)
	{
		double tmp = abs(v[n] - u[n]);
		d_uv += tmp*tmp;
	}

	return d_uv;
}
//****************End of || v - u||^2****************

vec_d utils::element_w_prod(const vec_d& v, const vec_d& u)
{
	
	vec_d out;

	for (int i = 0; i < v.size(); ++i) out.push_back(v[i] * u[i]);

	return out;

}

vec_d utils::const_prod(const vec_d& v, double w)
{
	vec_d out;

	for (int i = 0; i < v.size(); i++) out.push_back(v[i] * w);

	return out;
}

// interval limits update
void   update_LHS_interval(vec2_d& Large_HC_limits, vec2_d& HC_limits)
{
	int vec_num = Large_HC_limits[0].size();

	for (int r = 0; r < vec_num; r++)
	{
		if (Large_HC_limits[0][r] > HC_limits[0][r])
		{
			Large_HC_limits[0][r] = HC_limits[0][r];
		}
		if (Large_HC_limits[1][r] < HC_limits[1][r])
		{
			Large_HC_limits[1][r] = HC_limits[1][r];
		}
	}

}
//****************End of updating the interval limits****************

// Finding the center of each cluster
void   update_center_cluster(vec_d& Center_Pi, vec_d& Pi, int pi_num)
{
	int vec_num = Center_Pi.size();

	for (int r = 0; r < vec_num; r++)
	{
		Center_Pi[r] = Center_Pi[r] * ((double)(pi_num) / (pi_num + 1)) +
			((double)1 / (pi_num + 1)) * Pi[r];
	}
	

}
//****************End of Finding the center of each cluster****************

/// Mean VEctor
double utils::meanVector(const std::vector<double>& vec)
{
	double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
	double mean = sum / vec.size();
	return mean;
}

//****************End of mean vector****************

/// Mean VEctor
double utils::meanVector(const std::vector<int>& vec)
{
	double sum = std::accumulate(vec.begin(), vec.end(), 0.0);
	double mean = sum / vec.size();
	return mean;
}

//****************End of mean vector****************

/// STD VEctor
double utils::stdVector(const std::vector<double>& vec)
{
	double stdes = std::sqrt(utils::variance_func(vec));
	return stdes;
}

double utils::stdVector(const vec2_d& v)
{
	double stdes = 0;
	
	for (int i = 0; i < v.size(); i++)
	{
		stdes += std::sqrt(utils::variance_func(v[i]));
	}

	stdes = stdes / v.size();

	return stdes;
}

//****************End of STD vector****************

/// STD VEctor
double utils::stdVector(const std::vector<int>& vec)
{
	double stdes = std::sqrt(utils::variance_func(vec));
	return stdes;
}

//****************End of STD vector****************

/// Square VEctor
double utils::square_vec_d(const std::vector<double>& v)
{
	double sumv = accumulate(v.begin(), v.end(), 0.0);
	double aver = sumv / (double)v.size();
	double sqr = 0;
	for (int i = 0; i < v.size(); i++) sqr += (v[i] - aver)*(v[i] - aver);
	return sqr/(v.size() - 1);
}

/// Square VEctor
double utils::square_vec_d(const vec2_d& v)
{
	double ssw = 0;
	
	for (int i = 0; i < v.size(); i++)
	{
		ssw += square_vec_d(v[i]);
	}

	return ssw / (v.size() - 1);
}

//****************End of Square vector****************

void utils::print_vector(const vec_i& v)
{
	for (std::vector<int>::const_iterator i = v.begin(); i != v.end(); ++i)
		std::cout << *i << ' ';
}

void utils::print_vector(const vec_d& v)
{
	for (std::vector<double>::const_iterator i = v.begin(); i != v.end(); ++i)
		std::cout << *i << ' ';
	std::cout << std::endl;
}
void utils::print_vector(const vec2_d& v)
{
	for (vec2_d::const_iterator i = v.begin(); i != v.end(); ++i)
	{
		for (vec_d::const_iterator j = i->begin(); j != i->end(); ++j)
		{
			std::cout << *j << ' ';
		}
		std::cout << std::endl;
	}
		
	std::cout << std::endl;
}
//*****************************************************

double utils::sum_square(const vec_d& v)
{
	double m = utils::meanVector(v);
	double ss = 0;
	for (int i = 0; i < v.size(); i++) ss += (v[i] - m)*(v[i] - m);
	return ss;
}

//******************************************************

vec_d utils::get_Centre(const vec2_d& v)
{
	vec_d empt;

	for (int i = 0; i < v.size(); i++)
	{
		double tmp = 0;
		for (int j = 0; j < v[i].size(); j++)
		{
			tmp += v[i][j];
		}
		empt.push_back(tmp / (double)v[i].size());
	}

	return empt;
}

//******************************************************

double utils::L1_dist(const vec_d& v, const vec_d& u)
{
	double ans = 0;

	for (int i = 0; i < v.size(); i++)
	{
		ans += abs(v[i] - u[i]);
	}

	return ans;
}

//********************************************************

vec_d utils::prod_components(const vec_d& v, const vec_d& u)
{
	int vsize = v.size();
	int usize = u.size();
	vec_d empt;

	if (vsize == usize)
	{
		for (int i = 0; i < vsize; i++)
		{
			empt.push_back(v[i] * u[i]);
		}
	}
	else
	{
		printf("Prod Components: sizes are not equal!");
	}

	return empt;
}

//********************************************************

int utils::max_diff(const vec_d& v, double& val)
{
	int idx;
	int max_val = -INFINITY;

	for (int i = 0; i < v.size(); i++)
	{
		if (max_val < abs(v[i] - val))
		{
			idx = i;
			max_val = abs(v[i] - val);
 		}
	}

	return idx;

}

//********************************************************

int utils::max_diff(const vec_d& v, const vec_d& u)
{
	int idx;
	int max_val = -INFINITY;

	if (v.size() == u.size())
	{
		for (int i = 0; i < v.size(); i++)
		{
			if (max_val < abs(v[i] - u[i]))
			{
				max_val = abs(v[i] - u[i]);
				idx = i;
			}
		}
	}
	else
	{
		printf("Prod Components: sizes are not equal!");
	}


	return idx;
}

//********************************************************
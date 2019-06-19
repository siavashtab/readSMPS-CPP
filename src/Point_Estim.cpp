
#include "Point_Estim.h"

Point_Estim::Point_Estim()
{
	/*printf("\n");
	printf("--------------------------------------------------------------------------------\n");
	printf("                                   Point Estimation                             \n");
	printf("--------------------------------------------------------------------------------\n");*/

	//initialize the data structures
	cluster_samples = new std::vector<Scen_cluster>;
	x = new std::vector<SOL_str>;
	estimate_c = new vec_d;
	scen_idx = new vec2_i;
	dual_pool = new vec2_d;
	dual_pool_num = new vec_i;
	dual_pool_clust = new vec_i;
	samples = new vec2_d;
	samples_clust = new vec_i;
	scen_estim = new vec_d;

}

Point_Estim::~Point_Estim()
{
	cluster_samples->clear();
	std::vector<Scen_cluster>().swap(*cluster_samples);
	samples->clear();
	vec2_d().swap(*samples);
	samples_clust->clear();
	vec_i().swap(*samples_clust);
	dual_pool->clear();
	vec2_d().swap(*dual_pool);
	dual_pool_num->clear();
	vec_i().swap(*dual_pool_num);
	scen_idx->clear();
	vec2_i().swap(*scen_idx);
	scen_estim->clear();
	vec_d().swap(*scen_estim);
	estimate_c->clear();
	vec_d().swap(*estimate_c);
	dual_pool_clust->clear();
	dual_pool_clust->shrink_to_fit();
	estimate = 0;
	scen_num = 0;
}

void Point_Estim::end()
{
	cluster_samples->clear();
	std::vector<Scen_cluster>().swap(*cluster_samples);
	samples->clear();
	vec2_d().swap(*samples);
	samples_clust->clear();
	vec_i().swap(*samples_clust);
	dual_pool->clear();
	vec2_d().swap(*dual_pool);
	dual_pool_num->clear();
	vec_i().swap(*dual_pool_num);
	scen_idx->clear();
	vec2_i().swap(*scen_idx);
	scen_estim->clear();
	vec_d().swap(*scen_estim);
	estimate_c->clear();
	vec_d().swap(*estimate_c);
	dual_pool_clust->clear();
	dual_pool_clust->shrink_to_fit();
	estimate = 0;
	scen_num = 0;
}

void Point_Estim::initialize(int sampnum, int rep)
{
	
	estimate = 0;

	samp_num = sampnum;

	seed = 8879657642464524 + rep;

}


//Monte Carlo Point Estimation 
void  Point_Estim::monte(ProbPrep* problem, Sampling* sampl_ing)
{

	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());
	rng.seed(seed);
	std::srand(seed);


	scen_num = 0;
	initial_samp = samp_num;
	scen_idx->resize(initial_samp);

	std::vector<IloNum> Cx(problem->stage_sub_prob[0].range_raw.getSize(), -99900999);

	vec_d Ll;
	vec_d Ul;

	for (int rr = 0; rr < problem->SPprobINFO.RVs.size(); rr++)
	{
		Ll.push_back(0.0);
		Ul.push_back(1.0);
	}

	vec_d pi;

	for (int s = 0; s < initial_samp; s++)
	{

		Scen_cluster empt_clust;
		Scen_struct  empt_scen;


		//initialize
		empt_scen.samp = new vec_d;
		empt_scen.z1 = new vec_d;
		empt_scen.beta = new vec_d;
		empt_clust.samples = new vec2_d;
		empt_clust.size_of_samps = new vec_i;
		empt_clust.SA_vals = new vec2_d;
		empt_clust.SA_z1 = new vec2_d;
		empt_clust.z1 = new vec2_d;
		empt_clust.lower = new vec_d;
		empt_clust.upper = new vec_d;
		empt_clust.scens = new std::vector<Scen_struct>;
		empt_clust.beta_estim = new vec_d;
		empt_clust.samples_extrm = new vec2_d;
		empt_clust.pi = new vec_d;
		empt_clust.pi_tot = new vec_d;

		empt_clust.active = true;
		empt_clust.size = 1;
		empt_clust.clust_size = 1;
		empt_scen.active = true;

		//set the right hand side:
		sampl_ing->set_random_rhs_PE(*empt_scen.samp, *empt_scen.z1, problem->stage_sub_prob[0].rhs, problem->SPprobINFO.RVs,
			rng, Ll, Ul);
		samples->push_back(*empt_scen.samp);
		samples_clust->push_back(s);
		empt_clust.samples->push_back(*empt_scen.samp);
		empt_clust.size_of_samps->push_back(1);
		if(true) create_z1(*empt_clust.z1, *empt_scen.z1);

		//set the right hand side in the problem struct:
		if (problem->stage_sub_prob[0].beta.size() == 0)
		{
			problem->stage_sub_prob[0].beta.resize(problem->stage_sub_prob[0].range_raw.getSize());
			for (int r = 0; r < problem->stage_sub_prob[0].beta.size(); r++)
				problem->stage_sub_prob[0].beta[r].resize(problem->master_prob.num_var);
		}
		problem->solver.set_rhs_var(problem->stage_sub_prob[0]
			, *x, problem->stage_sub_prob[0].rhs, Cx);

		problem->solver.Solve_Prob(problem->stage_sub_prob[0], false);

		//fill out r_Cx
		int rv = 0;
		vec_d  r_Cx;
		vec_d Cxx;
		for (int r = 0; r < problem->stage_sub_prob[0].rhs.size(); r++)
		{
			if (problem->stage_sub_prob[0].rhs[r].isRandom == true)
			{
				r_Cx.push_back(empt_scen.samp->at(rv) - problem->stage_sub_prob[0].Cx[r]);
				Cxx.push_back(problem->stage_sub_prob[0].Cx[r]);
				rv++;
			}
		}

		//sual of the stochastic constraints
		pi = problem->stage_sub_prob[0].dual_R;
		empt_clust.r_Cx = r_Cx;
		empt_clust.Cx = Cxx;

		//Sensitivity analysis ranges for the RHS
		if (false)
		{
			*empt_clust.SA_vals = problem->solver.getRHSSA_rand(problem->stage_sub_prob[0],
				*empt_clust.lower, *empt_clust.upper);
			*empt_clust.SA_z1 = sampl_ing->inRHSSA(problem->SPprobINFO.RVs, *empt_clust.SA_vals,
				*empt_clust.lower, *empt_clust.upper);
			empt_clust.z1 = empt_clust.SA_z1;
		}
		

		//Add the discovered dual to the dual pool
		int pos;
		empt_clust.new_pi = !in2Vector(*dual_pool, pi, pos);
		if (empt_clust.new_pi == true || s == 0)
		{
			dual_pool->push_back(pi);
			dual_pool_num->push_back(1);
			empt_clust.pos = dual_pool_num->size() - 1;
			dual_pool_clust->push_back(s);
		}
		else
		{
			dual_pool_num->at(pos)++;
			empt_clust.pos = pos;
			if (false)
			{
				std::cout << "bef  " << std::endl;
				print_vectord(*cluster_samples->at(dual_pool_clust->at(pos)).z1);
			}

			if(true) update_z1(*cluster_samples->at(dual_pool_clust->at(pos)).z1, *empt_scen.z1);

			if (false)
			{
				std::cout << "af  " << std::endl;
				print_vectord(*cluster_samples->at(dual_pool_clust->at(pos)).z1);
				int ss = getchar();
			}
		}


		//scenario estimates:
		empt_scen.est = problem->stage_sub_prob[0].zstar;
		scen_estim->push_back(empt_scen.est);
		empt_clust.scen_idx.push_back(s);
		empt_clust.scens->push_back(empt_scen);


		//Cut Coeffs
		empt_scen.alpha = std::inner_product(problem->stage_sub_prob[0].dual_tot.begin(),
			problem->stage_sub_prob[0].dual_tot.end(), problem->stage_sub_prob[0].r_w.begin(), 0.0);
		empt_scen.beta->resize(problem->master_prob.num_var);
		for (int r = 0; r < problem->stage_sub_prob[0].beta.size(); r++)
		{
			for (int v = 0; v < problem->master_prob.num_var; v++)
			{
				empt_scen.beta->at(v) += problem->stage_sub_prob[0].beta[r][v] * problem->stage_sub_prob[0].dual_tot[r];
			}
		}
		
		
		//Estimation
		empt_clust.estim = empt_scen.est;
		empt_clust.alpha_estim = empt_scen.alpha;
		empt_clust.beta_estim = empt_scen.beta;
		empt_clust.estim_R = std::inner_product(pi.begin(), pi.end(), problem->stage_sub_prob[0].rho.begin(), 0.0);
		empt_clust.estim_D = empt_clust.estim - empt_clust.estim_R;
		if (false)
		{
			empt_clust.estim_lower = empt_clust.estim_D + std::inner_product(pi.begin(), pi.end(), empt_clust.lower->begin(), 0.0);
			empt_clust.estim_upper = empt_clust.estim_D + std::inner_product(pi.begin(), pi.end(), empt_clust.upper->begin(), 0.0);
			empt_clust.estim_extrm.push_back(empt_clust.estim_lower);
			empt_clust.estim_extrm.push_back(empt_clust.estim_upper);
			empt_clust.samples_extrm->push_back(*empt_clust.lower);
			empt_clust.samples_extrm->push_back(*empt_clust.upper);
		}
		empt_clust.SSW = 0;
		empt_clust.SSW = square_vec_d(empt_clust.estim_extrm);
		
		//Check the solution and estimates
		if (false)
		{
			problem->solver.Error_Report(problem->stage_sub_prob[0], "ERROR: Point Estim: estimR > estim");
			printf("r - Cx:  \n");
			print_vectord(r_Cx);
			printf("Pi:  \n");
			print_vectord(pi);
			printf("x:  \n");
			for (int vv = 0; vv < x->size(); vv++) std::cout << x->at(vv).col << " " << x->at(vv).value << std::endl;
			printf("Estim R:  %0.04f  - Estim:  %0.04f \n", empt_clust.estim_R, empt_clust.estim);
			int ss = getchar();
		}

	
		//Pi is the like a tag for each cluster
		*empt_clust.pi = pi;
		double rmax_ = maxVector(pi);
		empt_clust.maxr_ = rmax_;
		*empt_clust.pi_tot = problem->stage_sub_prob[0].dual_tot;
		
		if (print_SSW_info == 1)
		{
			std::cout << "--------------------------------" << std::endl;
			std::cout << "cluster:  " << s << std::endl;
			print_vectord(empt_clust.samples->at(0));
			std::cout << std::endl;
			std::cout << "dual: " << std::endl;
			print_vectord(pi);
			std::cout << std::endl << std::endl;
			print_vectord(*empt_clust.lower);
			std::cout << std::endl;
			print_vectord(*empt_clust.upper);
			std::cout << std::endl << std::endl;
			std::cout << "lower : " << empt_clust.estim_lower << "  -  upper : " <<
				empt_clust.estim_upper << std::endl;
			std::cout << "SSW:   " << empt_clust.SSW << std::endl;
			std::cout << "new pi:   " << std::boolalpha << empt_clust.new_pi << std::endl;
			if (empt_clust.new_pi == false)
			{
				//cout << dual_pool.size() << " " << empt_clust.pos << endl;
				print_vectord(dual_pool->at(empt_clust.pos));
				std::cout << std::endl;
				std::cout << dual_pool_num->at(empt_clust.pos) << std::endl;
			}
			std::cout << "--------------------------------" << std::endl;
		}

		//Add to the pool of clusters
		cluster_samples->push_back(empt_clust);

		//Clear the created structs
		(*scen_idx)[s].push_back(s);
		pi.clear();
		
		//one sample is added
		scen_num++;

	}

	if (print_SSW_info == 1)
	{
		int ss = getchar();
	}

	//Update the estimators
	update_estimate();


	Ll.clear();
	Ll.shrink_to_fit();
	Ul.clear();
	Ul.shrink_to_fit();
	pi.clear();
	pi.shrink_to_fit();
	Cx.clear();
	Cx.shrink_to_fit();

	bool sampprint = false;
	if (sampprint == true)
	{
		for (int cl = 0; cl < cluster_samples->size(); cl++)
		{
			printf("before cl %d  -  size: %d  -  estim:  &%0.02f \n", 
				cl, cluster_samples->at(cl).size, cluster_samples->at(cl).estim);
			print_vectord(*cluster_samples->at(cl).pi);
			print_vectord(*cluster_samples->at(cl).scens->at(0).samp);
			printf("\n\n");
		}
	}

}
//**************************************************************************

//LHS Point Estimation
void  Point_Estim::lhs(ProbPrep* problem, Sampling* sampl_ing)
{

	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());
	rng.seed(seed);
	std::srand(seed);

	scen_num = 0;
	initial_samp = samp_num;
	scen_idx->resize(initial_samp);


	std::vector<IloNum> Cx(problem->stage_sub_prob[0].range_raw.getSize(), -99900999);

	double interval = 1 / (double)initial_samp;

	std::vector<int> perm;
	vec2_i  perms;
	for (int h = 0; h < initial_samp; ++h) perm.push_back(h + 1); // 1 2 3 4 5 6 7 8 9

																  //Generate Permutations for each random variable
	for (int rr = 0; rr < problem->SPprobINFO.RVs.size(); rr++)
	{
		std::random_shuffle(perm.begin(), perm.end());
		perms.push_back(perm);
	}
	perm.clear();

	vec_d pi;

	for (int s = 0; s < initial_samp; s++)
	{

		// Creating the 0-1 intervals for LHS *************************
		vec_d Ll;
		vec_d Ul;
		vec_i permi;

		//Find the 0-1 interval associated with each permutation
		for (int rr = 0; rr < problem->SPprobINFO.RVs.size(); rr++)
		{
			Ll.push_back((perms[rr][s] - 1) * interval);
			Ul.push_back(perms[rr][s] * interval);
			permi.push_back(perms[rr][s]);
		}
		//**************************************************************

		Scen_cluster empt_clust;
		Scen_struct  empt_scen;

		empt_clust.active = true;
		empt_clust.size = 1;
		empt_clust.clust_size = 1;
		empt_scen.active = true;

		//initialize
		empt_scen.samp = new vec_d;
		empt_scen.z1 = new vec_d;
		empt_scen.beta = new vec_d;
		empt_clust.samples = new vec2_d;
		empt_clust.size_of_samps = new vec_i;
		empt_clust.SA_vals = new vec2_d;
		empt_clust.SA_z1 = new vec2_d;
		empt_clust.z1 = new vec2_d;
		empt_clust.lower = new vec_d;
		empt_clust.upper = new vec_d;
		empt_clust.scens = new std::vector<Scen_struct>;
		empt_clust.beta_estim = new vec_d;
		empt_clust.samples_extrm = new vec2_d;
		empt_clust.pi = new vec_d;
		empt_clust.pi_tot = new vec_d;

		//set the right hand side:
		sampl_ing->set_random_rhs_PE(*empt_scen.samp, *empt_scen.z1, problem->stage_sub_prob[0].rhs, problem->SPprobINFO.RVs,
			rng, Ll, Ul);
		samples->push_back(*empt_scen.samp);
		samples_clust->push_back(s);
		empt_clust.samples->push_back(*empt_scen.samp);
		empt_clust.size_of_samps->push_back(1);
		create_z1(*empt_clust.z1, *empt_scen.z1);

		//set the right hand side in the problem struct:
		if (problem->stage_sub_prob[0].beta.size() == 0)
		{
			problem->stage_sub_prob[0].beta.resize(problem->stage_sub_prob[0].range_raw.getSize());
			for (int r = 0; r < problem->stage_sub_prob[0].beta.size(); r++)
				problem->stage_sub_prob[0].beta[r].resize(problem->master_prob.num_var);
		}
		problem->solver.set_rhs_var(problem->stage_sub_prob[0]
			, *x, problem->stage_sub_prob[0].rhs, Cx);

		problem->solver.Solve_Prob(problem->stage_sub_prob[0], false);

		//fill out r_Cx
		int rv = 0;
		vec_d  r_Cx;
		vec_d  Cxx;
		for (int r = 0; r < problem->stage_sub_prob[0].rhs.size(); r++)
		{
			if (problem->stage_sub_prob[0].rhs[r].isRandom == true)
			{
				r_Cx.push_back(empt_scen.samp->at(rv) - problem->stage_sub_prob[0].Cx[r]);
				rv++;
			}
			Cxx.push_back(problem->stage_sub_prob[0].Cx[r]);
		}

		//sual of the stochastic constraints
		pi = problem->stage_sub_prob[0].dual_R;
		empt_clust.r_Cx = r_Cx;
		empt_clust.Cx = Cxx;

		//Sensitivity analysis ranges for the RHS
		//*empt_clust.SA_vals = problem->solver.getRHSSA_rand(problem->stage_sub_prob[0],
		//	*empt_clust.lower, *empt_clust.upper);
		//*empt_clust.SA_z1 = sampl_ing->inRHSSA(problem->SPprobINFO.RVs, *empt_clust.SA_vals,
		//	*empt_clust.lower, *empt_clust.upper);
		//empt_clust.z1 = empt_clust.SA_z1;

		//Add the discovered dual to the dual pool
		int pos;
		empt_clust.new_pi = !in2Vector(*dual_pool, pi, pos);
		if (empt_clust.new_pi == true || s == 0)
		{
			dual_pool->push_back(pi);
			dual_pool_num->push_back(1);
			empt_clust.pos = dual_pool_num->size() - 1;
			dual_pool_clust->push_back(s);
		}
		else
		{
			dual_pool_num->at(pos)++;
			empt_clust.pos = pos;
			if (false)
			{
				std::cout << "bef  " << std::endl;
				print_vectord(*cluster_samples->at(dual_pool_clust->at(pos)).z1);
			}

			if (true) update_z1(*cluster_samples->at(dual_pool_clust->at(pos)).z1, *empt_scen.z1);

			if (false)
			{
				std::cout << "af  " << std::endl;
				print_vectord(*cluster_samples->at(dual_pool_clust->at(pos)).z1);
				int ss = getchar();
			}
		}


		//scenario estimates:
		empt_scen.est = problem->stage_sub_prob[0].zstar;
		empt_clust.scen_idx.push_back(s);


		//Cut Coeffs
		empt_scen.alpha = inner_product(problem->stage_sub_prob[0].dual_tot.begin(),
			problem->stage_sub_prob[0].dual_tot.end(), problem->stage_sub_prob[0].r_w.begin(), 0.0);
		empt_scen.beta->resize(problem->master_prob.num_var);

		empt_clust.scens->push_back(empt_scen);

		for (int r = 0; r < problem->stage_sub_prob[0].beta.size(); r++)
		{
			for (int v = 0; v < problem->master_prob.num_var; v++)
			{
				empt_scen.beta->at(v) += problem->stage_sub_prob[0].beta[r][v] * problem->stage_sub_prob[0].dual_tot[r];
			}
		}

		//Estimation
		empt_clust.estim = empt_scen.est;
		scen_estim->push_back(empt_scen.est);
		empt_clust.alpha_estim = empt_scen.alpha;
		empt_clust.beta_estim = empt_scen.beta;
		empt_clust.estim_R = inner_product(pi.begin(), pi.end(), problem->stage_sub_prob[0].rho.begin(), 0.0);
		empt_clust.estim_D = empt_clust.estim - empt_clust.estim_R;
		//empt_clust.estim_lower = empt_clust.estim_D + inner_product(pi.begin(), pi.end(), empt_clust.lower->begin(), 0.0);
		//empt_clust.estim_upper = empt_clust.estim_D + inner_product(pi.begin(), pi.end(), empt_clust.upper->begin(), 0.0);
		//empt_clust.estim_extrm.push_back(empt_clust.estim_lower);
		//empt_clust.estim_extrm.push_back(empt_clust.estim_upper);
		//empt_clust.samples_extrm->push_back(*empt_clust.lower);
		//empt_clust.samples_extrm->push_back(*empt_clust.upper);
		empt_clust.SSW = 0;
		//empt_clust.SSW = square_vec_d(empt_clust.estim_extrm);

		//Check the solution and estimates
		if (false)
		{
			printf("r(w): scen: %d \n", s);
			print_vectord(*empt_scen.samp);
			problem->solver.Error_Report(problem->stage_sub_prob[0], "ERROR: Point Estim: estimR > estim");
			printf("r - Cx:  \n");
			print_vectord(r_Cx);
			printf("Pi:  \n");
			print_vectord(pi);
			printf("x:  \n");
			for (int vv = 0; vv < x->size(); vv++) std::cout << x->at(vv).col << " " << x->at(vv).value << std::endl;
			printf("Estim R:  %0.04f  - Estim:  %0.04f \n", empt_clust.estim_R, empt_clust.estim);
			int ss = getchar();
		}


		//Pi is the like a tag for each cluster
		*empt_clust.pi = pi;
		double rmax_ = maxVector(pi);
		empt_clust.maxr_ = rmax_;
		*empt_clust.pi_tot = problem->stage_sub_prob[0].dual_tot;

		if (print_SSW_info == 1)
		{
			std::cout << "--------------------------------" << std::endl;
			std::cout << "cluster:  " << s << std::endl;
			print_vectord(empt_clust.samples->at(0));
			std::cout << std::endl;
			std::cout << "dual: " << std::endl;
			print_vectord(pi);
			std::cout << std::endl << std::endl;
			print_vectord(*empt_clust.lower);
			std::cout << std::endl;
			print_vectord(*empt_clust.upper);
			std::cout << std::endl << std::endl;
			std::cout << "lower : " << empt_clust.estim_lower << "  -  upper : " <<
				empt_clust.estim_upper << std::endl;
			std::cout << "SSW:   " << empt_clust.SSW << std::endl;
			std::cout << "new pi:   " << std::boolalpha << empt_clust.new_pi << std::endl;
			if (empt_clust.new_pi == false)
			{
				//cout << dual_pool.size() << " " << empt_clust.pos << endl;
				print_vectord(dual_pool->at(empt_clust.pos));
				std::cout << std::endl;
				std::cout << dual_pool_num->at(empt_clust.pos) << std::endl;
			}
			std::cout << "--------------------------------" << std::endl;
		}

		//Add to the pool of clusters
		cluster_samples->push_back(empt_clust);

		//Clear the created structs
		(*scen_idx)[s].push_back(s);

		//one sample is added
		scen_num++;

		Ll.clear();
		Ll.shrink_to_fit();
		Ul.clear();
		Ul.shrink_to_fit();
		permi.clear();
		permi.shrink_to_fit();

	}

	update_estimate();


	//cout << mean_MSW << " " << estimate_std << " " << 1 - (mean_MSW / pow(estimate_std, 2)) << endl;

	if (print_SSW_info == 1)
	{
		int ss = getchar();
	}

	perm.clear();
	perm.shrink_to_fit();
	perms.clear();
	perms.shrink_to_fit();
	pi.clear();
	pi.shrink_to_fit();
	Cx.clear();
	Cx.shrink_to_fit();
	

}
//**************************************************************************


//ACS Point Estimation
void  Point_Estim::acs(ProbPrep* problem, Sampling* sampl_ing)
{
	
	dual_num = 0;

	if (scen_or_clust == 0)
	{
		if (scen_type == 0)
		{
			monte(problem, sampl_ing);
		}
		else
		{
			lhs(problem, sampl_ing);
		}
	}
	else
	{
		seq(problem, sampl_ing);
	}
	
	
	isACS = true;
	int ult_size = samp_num + max_n_c;
	double R_a_prev = R_a;

	int iter = 0;
	//If it is benificial to do ACS do the following
	if (true)
	{

		double est_copy = estimate;

		while (R_a > 0 && iter < 5)
		{

			R_a_prev = R_a;

			for (int cl = 0; cl < cluster_samples->size(); cl++)
			{

				int size_ = cluster_samples->at(cl).size;

				if (n_c_tot < 1) n_c_tot = 1;
				int num = n_c_tot * 4;


				if (cluster_samples->at(cl).active == true&&
					cluster_samples->at(cl).new_pi == true)
				{
					if (false)
					{
						printf("iter: %d - clust: %d  -  how many duals?: %d  \n", iter, cl, dual_num);
					}
					if (print_point_Estim_samples == 1)
					{
						printf("before: cluster %d , \n", cl);

						for (int i = 0; i < cluster_samples->at(cl).samples->size(); i++)
						{
							print_vectord(cluster_samples->at(cl).samples->at(i));
							printf("  estim: %0.06f \n", cluster_samples->at(cl).scens->at(i).est);
						}

						printf("before: dual , \n");
						print_vectord(*cluster_samples->at(cl).pi);
						printf("before: SSW mean: %0.06f , \n", mean_MSW);

						printf("before: estimate of %d: %0.06f , SSW: %0.06f , tot estimate: %0.06f , \n"
							, cl, cluster_samples->at(cl).estim, cluster_samples->at(cl).SSW, estimate);

						printf(" -------------------------------- \n \n", cl);
					}

					try
					{
						//create_neighbor_None(problem, sampl_ing, cl, num);
						//create_neighbor_None_rmax(problem, sampl_ing, cl, num);
						//create_neighbor_all(problem, cl);
						create_neighbor(problem, sampl_ing, cl, num);
					}
					catch (...)
					{
						std::cerr << "ERROR: Point Estim: ACS: Create Neighbor" << std::endl;
						system("pause");
					}
					
					if (scen_or_clust == 0)
					{
						update_estimate();
					}
					else
					{
						update_estimate_seq();
					}
					
					

					if (print_point_Estim_samples == 1)
					{
						printf("after: cluster %d , \n", cl);

						for (int i = 0; i < cluster_samples->at(cl).samples->size(); i++)
						{
							print_vectord(cluster_samples->at(cl).samples->at(i));
							printf("  estim: %0.06f \n", cluster_samples->at(cl).scens->at(i).est);
						}

						printf("after: estimate of %d: %0.06f , SSW: %0.06f , tot estimate: %0.06f , \n"
							, cl, cluster_samples->at(cl).estim, cluster_samples->at(cl).SSW, estimate);

						printf(" ------------------------------------------------------------- \n ", cl);
						printf(" ------------------------------------------------------------- \n \n ", cl);

						int ss = getchar();
					}



				}
			}

			iter++;

			if (R_a_prev == R_a) break;
		}


	}

	bool sampprint = false;
	if (sampprint == true)
	{
		for (int cl = 0; cl < cluster_samples->size(); cl++)
		{
			printf("after cl %d  -  size: %d  -  estim:  &%0.02f \n",
				cl, cluster_samples->at(cl).size, cluster_samples->at(cl).estim);
			for (int i = 0; i < cluster_samples->at(cl).scens->size(); i++)
				print_vectord(*cluster_samples->at(cl).scens->at(i).samp);
			printf("\n\n");
		}
		int ss = getchar();
	}

	//system("pause");

}
//**************************************************************************

//**************************************************************************


//**************************************************************************

//**************************************************************************
//Point estimate across all clusters
void Point_Estim::update_estimate()
{
	estimate_c->clear();
	estimate = 0;
	mean_MSW = 0;
	int newPis = 0;
	bool printing = false;

	for (int i = 0; i < cluster_samples->size(); i++)
	{
		estimate += (1 / (double)cluster_samples->size()) * cluster_samples->at(i).estim;

		
		estimate_c->push_back(cluster_samples->at(i).estim);
		double SSW = 0;
		for (int j = 0; j < cluster_samples->at(i).scens->size(); j++)
			SSW += (cluster_samples->at(i).scens->at(j).est - cluster_samples->at(i).estim)*
			                              (cluster_samples->at(i).scens->at(j).est - cluster_samples->at(i).estim);
		if (cluster_samples->at(i).SSW_it.size() > 1)
			if (cluster_samples->at(i).SSW_it[cluster_samples->at(i).SSW_it.size() - 1] - SSW > 0.01*SSW)
				cluster_samples->at(i).active = false;
		cluster_samples->at(i).SSW_it.push_back(SSW);
		cluster_samples->at(i).SSW = SSW;

		if (printing == true)
		{
			printf("cluster %d  ,  %d , %0.4f , %0.02f \n", i, cluster_samples->at(i).size, cluster_samples->at(i).estim,
				SSW);
			std::cout << "clust active: " << cluster_samples->at(i).active << std::endl;
			std::cout << estimate << std::endl;
		}
	}


	estimate_std = std_vec_d(*estimate_c);

	n_c_tot = 0;
	mean_MSW = 0;
	for (int i = 0; i < cluster_samples->size(); i++)
	{
		if (cluster_samples->at(i).new_pi == true)
		{
			newPis++;
			mean_MSW += cluster_samples->at(i).SSW ;
		}		
	}

	mean_MSW = mean_MSW / (double)newPis;

	MSW_it.push_back(mean_MSW);
	n_c_tot = ceil(sqrt((1 - R_a)*(1 - R_a) / (R_a * R_a)));
	R_a = 1 - (mean_MSW / pow(estimate_std,2));

	if (printing == true)
	{
		printf("MSW = %0.02f - n_c = %d  - R_a = %0.02f \n", mean_MSW, n_c_tot, R_a);
		int ss = getchar();
	}

	dual_num = dual_pool->size();

}
//**************************************************************************

//**************************************************************************
//Point estimate across all clusters
void Point_Estim::update_estimate_seq()
{
	
	estimate_c->clear();
	estimate = 0;
	mean_MSW = 0;
	int counter = 0;
	double sum_est = 0;
	bool printing = false;
	int newPis = 0;

	for (int i = 0; i < cluster_samples->size(); i++)
	{
		if (cluster_samples->at(i).new_pi == true)
		{
			//counter += dual_pool_num->at(cluster_samples->at(i).pos);
			//counter += cluster_samples->at(i).clust_size;
			counter += cluster_samples->at(i).clust_size;
			sum_est += cluster_samples->at(i).estim * cluster_samples->at(i).clust_size;
			if (printing == true)std::cout << cluster_samples->at(i).clust_size << "  ";
			if (printing == true)std::cout << cluster_samples->at(i).estim << std::endl;
			estimate_c->push_back(cluster_samples->at(i).estim);
		}
	}

	estimate = sum_est / (double)counter;
	estimate_std = std_vec_d(*estimate_c);

	n_c_tot = 0;
	mean_MSW = 0;
	for (int i = 0; i < cluster_samples->size(); i++)
	{
		if (cluster_samples->at(i).new_pi == true)
		{
			newPis++;
			mean_MSW += cluster_samples->at(i).SSW;
		}
	}

	mean_MSW = mean_MSW / (double)newPis;

	MSW_it.push_back(mean_MSW);
	n_c_tot = ceil(sqrt((1 - R_a)*(1 - R_a) / (R_a * R_a)));
	R_a = 1 - (mean_MSW / pow(estimate_std, 2));

	if (printing == true)
	{
		printf("MSW = %0.02f - n_c = %d  - R_a = %0.02f \n", mean_MSW, n_c_tot, R_a);
		int ss = getchar();
	}

	dual_num = dual_pool->size();

}
//**************************************************************************

//**************************************************************************
// create neighbor samples
void Point_Estim::create_neighbor(ProbPrep* problem, Sampling* sampl_ing, int clust, int num)
{

	std::vector<IloNum> Cx(problem->stage_sub_prob[0].num_rng, -99900999);

	// Secondary set of samples:
	//std::vector<Scen_struct> poten_scens = monte2(problem, sampl_ing, *cluster_samples->at(clust).z1, num);
	std::vector<Scen_struct> poten_scens = lhs2(problem, sampl_ing, *cluster_samples->at(clust).z1, num);
	
	vec2_d zero1 = *cluster_samples->at(clust).SA_z1;
	bool printing = false;
	
	for (int s = 0; s < poten_scens.size(); s++)
	{
		if (true)
		{
			int pos;

			if (in2Vector(*samples, *poten_scens[s].samp, pos) == false)
			{
				
				bool isNew = true;

				////Check if the generated sample is new____________________________________________
				for (int sce = 0; sce < cluster_samples->at(clust).scens->size(); sce++)
				{
					if (*poten_scens[s].samp == *cluster_samples->at(clust).scens->at(sce).samp)
					{
						cluster_samples->at(clust).size++;

						//cut coeffs
						cluster_samples->at(clust).alpha_estim =
							((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
							cluster_samples->at(clust).alpha_estim +
							(1 / (double)cluster_samples->at(clust).size) *
							cluster_samples->at(clust).scens->at(sce).alpha;

						//Estimation
						cluster_samples->at(clust).estim =
							((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
							cluster_samples->at(clust).estim +
							(1 / (double)cluster_samples->at(clust).size) *
							cluster_samples->at(clust).scens->at(sce).est;

						isNew = false;
						break;
					}
				}
				////_____________________________________________________________________________________

				if (isNew == true)
				{
					//Put the sample in rhs struct
					int rv = 0;
					for (int r = 0; r < problem->stage_sub_prob[0].rhs.size(); r++)
					{
						if (problem->stage_sub_prob[0].rhs[r].isRandom == true)
						{
							problem->stage_sub_prob[0].rhs[r].value = (IloNum)poten_scens[s].samp->at(rv);
							rv++;
						}
					}

					problem->solver.set_rhs_var(problem->stage_sub_prob[0]
						, *x, problem->stage_sub_prob[0].rhs, Cx);

					problem->solver.Solve_Prob(problem->stage_sub_prob[0], false);

					//Check the new generated samples
					if (printing == true)
					{
						printf("cluster:  %d   \n", clust);
						printf("current sample   \n");
						print_vectord(cluster_samples->at(clust).samples->at(0));
						printf("cluster size: %d   \n", cluster_samples->at(clust).size);
						std::cout << "new pi:  " << cluster_samples->at(clust).new_pi << std::endl;
						printf("pi val:   \n");
						print_vectord(*cluster_samples->at(clust).pi);
						std::cout << std::endl << std::endl;
						for (int r = 0; r < poten_scens[s].samp->size(); r++)
						{
							printf("lower : %0.02f - sample: %0.02f - upper: %0.02f  \n",
								cluster_samples->at(clust).SA_vals->at(r)[0],
								poten_scens[s].samp->at(r),
								cluster_samples->at(clust).SA_vals->at(r)[1]);
						}
						printf("estimates   \n");
						printf("clust: %d  -  scen estim: %0.04f   -   cluster estim: %0.04f  \n", clust,
							poten_scens[s].est, cluster_samples->at(clust).estim);
						printf("estimate   \n");
						double est = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
							cluster_samples->at(clust).estim +
							(1 / (double)cluster_samples->at(clust).size) * poten_scens[s].est;
						std::cout << est << std::endl;
						printf("check its dual:   \n");
						print_vectord(problem->stage_sub_prob[0].dual_R);
						printf("---------------------------------------------------------------------- \n\n");
						int ss = getchar();
					}

					for (int cl = 0; cl < cluster_samples->size(); cl++)
					{
						if (cluster_samples->at(cl).new_pi == true)
						{
							if (problem->stage_sub_prob[0].dual_R == *cluster_samples->at(cl).pi)
							{
								clust = cl;
								//samples->push_back(*poten_scens[s].samp);
								samples_clust->push_back(cl);
								cluster_samples->at(cl).size++;
								cluster_samples->at(cl).samples->push_back(*poten_scens[s].samp);
								cluster_samples->at(cl).size_of_samps->push_back(1);
								update_z1(*cluster_samples->at(cl).z1, *poten_scens[s].z1);

								poten_scens[s].est = problem->stage_sub_prob[0].zstar;
								scen_estim->push_back(poten_scens[s].est);
								cluster_samples->at(cl).scen_idx.push_back(scen_num);
								cluster_samples->at(cl).scens->push_back(poten_scens[s]);

								//Cut Coeffs
								poten_scens[s].alpha = inner_product(cluster_samples->at(cl).pi_tot->begin(),
									cluster_samples->at(cl).pi_tot->end(), problem->stage_sub_prob[0].r_w.begin(), 0.0);


								cluster_samples->at(cl).alpha_estim =
									((cluster_samples->at(cl).size - 1) / (double)cluster_samples->at(cl).size) *
									cluster_samples->at(cl).alpha_estim +
									(1 / (double)cluster_samples->at(cl).size) * poten_scens[s].alpha;

								//Estimation
								cluster_samples->at(cl).estim = ((cluster_samples->at(cl).size - 1) / (double)cluster_samples->at(cl).size) *
									cluster_samples->at(cl).estim +
									(1 / (double)cluster_samples->at(cl).size) * poten_scens[s].est;

								scen_idx->at(cl).push_back(scen_num);

								double SSW = 0;
								for (int sc = 0; sc < cluster_samples->at(cl).scens->size(); sc++)
								{
									SSW += (cluster_samples->at(cl).scens->at(sc).est - cluster_samples->at(cl).estim)*
										(cluster_samples->at(cl).scens->at(sc).est - cluster_samples->at(cl).estim);
								}
								if (SSW <= cluster_samples->at(cl).SSW) cluster_samples->at(cl).active = false;
								cluster_samples->at(cl).SSW = SSW;
								scen_num++;

								break;
							}//pi is equal
						}
					}

				}



			}		
		}
	}

	poten_scens.clear();
	poten_scens.shrink_to_fit();
	Cx.clear();
	Cx.shrink_to_fit();

}
//**************************************************************************

//**************************************************************************
//Check if the obtained interval is inside the SA range
bool Point_Estim::is_in_SArange(vec2_d& SA, vec_d& z1)
{
	bool output = true;

	if (z1.size() != SA.size())
	{
		printf("Error Point_Estim: is in Range: z1 size: %d , SA size: %d \n", z1.size(), SA.size());
		system("pause");
	}
		

	for (int i = 0; i < z1.size(); i++)
	{
		//cout << SA[i][0] << " " << z1[i] << " " << SA[i][1] << endl;
		if (z1[i] > SA[i][1] || z1[i] < SA[i][0])
		{
			output = false;
			break;
		}
	}
	
	//system("pause");

	return output;

}
//**************************************************************************

//**************************************************************************
//Secondary LHS 
std::vector<Scen_struct>  Point_Estim::lhs2(ProbPrep* problem, Sampling* sampl_ing, vec2_d z1, int size)
{
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());
	rng.seed(seed);
	std::srand(seed);

	std::vector<Scen_struct> output;

	std::vector<int> perm;
	vec2_i  perms;
	for (int h = 0; h < size; ++h) perm.push_back(h + 1); // 1 2 3 4 5 6 7 8 9

																  //Generate Permutations for each random variable
	for (int rr = 0; rr < problem->SPprobINFO.RVs.size(); rr++)
	{
		std::random_shuffle(perm.begin(), perm.end());
		perms.push_back(perm);
	}
	perm.clear();


	for (int s = 0; s < size; s++)
	{

		Scen_struct  empt_scen;
		empt_scen.samp = new vec_d;
		empt_scen.z1   = new vec_d;
		empt_scen.beta = new vec_d;

		vec_d Ll;
		vec_d Ul;
		vec_i permi;

		//Find the 0-1 interval associated with each permutation
		for (int rr = 0; rr < problem->SPprobINFO.RVs.size(); rr++)
		{
			//printf("scen: %d , rv: , %d , %0.06f , %0.06f \n", s, rr, z1[rr][0], z1[rr][1]);

			double interval = abs(((double)(z1[rr][1] - z1[rr][0])) / size);
			Ll.push_back(std::min(1.0 - ACS_cov,z1[rr][0] + (perms[rr][s] - 1) * interval));
			Ul.push_back(std::min(1.0,z1[rr][0] + perms[rr][s] * interval));
			//Ll.push_back(min(1.0 - ACS_cov, 0 + (perms[rr][s] - 1) * interval));
			//Ul.push_back(min(1.0, 0 + perms[rr][s] * interval));
		}


		empt_scen.active = true;

		//set the right hand side:
		sampl_ing->set_random_rhs_PE(*empt_scen.samp, *empt_scen.z1, problem->stage_sub_prob[0].rhs, problem->SPprobINFO.RVs,
			rng, Ll, Ul);
		
		output.push_back(empt_scen);


		Ll.clear();
		Ul.clear();


	}


	perm.clear();
	perms.clear();

	return output;

}
//**************************************************************************

//**************************************************************************
//Secondary SRS
std::vector<Scen_struct>  Point_Estim::monte2(ProbPrep* problem, Sampling* sampl_ing, vec2_d z1, int size)
{
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());
	rng.seed(seed);
	std::srand(seed);

	std::vector<Scen_struct> output;
	vec_d Ll;
	vec_d Ul;

	//Find the 0-1 interval associated with each permutation
	for (int rr = 0; rr < problem->SPprobINFO.RVs.size(); rr++)
	{
		double interval = abs(((double)(z1[rr][1] - z1[rr][0])) / size);
		Ll.push_back(z1[rr][0]);
		Ul.push_back(z1[rr][1]);
	}


	for (int s = 0; s < size; s++)
	{

		Scen_struct  empt_scen;
		empt_scen.samp = new vec_d;
		empt_scen.z1 = new vec_d;
		empt_scen.beta = new vec_d;
		empt_scen.active = true;

		//set the right hand side:
		sampl_ing->set_random_rhs_PE(*empt_scen.samp, *empt_scen.z1, problem->stage_sub_prob[0].rhs, problem->SPprobINFO.RVs,
			rng, Ll, Ul);

		output.push_back(empt_scen);


	}

	Ll.clear();
	Ll.shrink_to_fit();
	Ul.clear();
	Ul.shrink_to_fit();


	return output;

}
//**************************************************************************

//Sequential sampling:
//Monte Carlo Point Estimation 
void  Point_Estim::seq(ProbPrep* problem, Sampling* sampl_ing)
{
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());
	rng.seed(seed);
	std::srand(seed);


	int min_iter = 100;
	int max_iter = 1000;
	int iter = 0;
	tot_size = 0;

	scen_num = 0;
	initial_samp = samp_num;
	scen_idx->resize(initial_samp);

	std::vector<IloNum> Cx(problem->stage_sub_prob[0].range_raw.getSize(), -99900999);

	vec_d Ll;
	vec_d Ul;

	for (int rr = 0; rr < problem->SPprobINFO.RVs.size(); rr++)
	{
		Ll.push_back(0.0);
		Ul.push_back(1.0);
	}

	vec_d pi;
	mean_MSW = 0;
	
	while(iter  < samp_num)
	{

		Scen_cluster empt_clust;
		Scen_struct  empt_scen;


		empt_clust.active = true;
		empt_clust.size = 1;
		empt_scen.active = true;

		//initialize
		empt_scen.samp = new vec_d;
		empt_scen.z1 = new vec_d;
		empt_scen.beta = new vec_d;
		empt_clust.samples = new vec2_d;
		empt_clust.size_of_samps = new vec_i;
		empt_clust.SA_vals = new vec2_d;
		empt_clust.SA_z1 = new vec2_d;
		empt_clust.z1 = new vec2_d;
		empt_clust.lower = new vec_d;
		empt_clust.upper = new vec_d;
		empt_clust.scens = new std::vector<Scen_struct>;
		empt_clust.beta_estim = new vec_d;
		empt_clust.samples_extrm = new vec2_d;
		empt_clust.pi = new vec_d;
		empt_clust.pi_tot = new vec_d;
		
		//set the right hand side:
		sampl_ing->set_random_rhs_PE(*empt_scen.samp, *empt_scen.z1, problem->stage_sub_prob[0].rhs, problem->SPprobINFO.RVs,
			rng, Ll, Ul);
		
		int samp_pos;
		//if (in2Vector(*samples, *empt_scen.samp, samp_pos) == false)
		if (true)
		{
			empt_clust.samples->push_back(*empt_scen.samp);
			empt_clust.size_of_samps->push_back(1);
			samples->push_back(*empt_scen.samp);

			//set the right hand side in the problem struct:
			if (problem->stage_sub_prob[0].beta.size() == 0)
			{
				problem->stage_sub_prob[0].beta.resize(problem->stage_sub_prob[0].range_raw.getSize());
				for (int r = 0; r < problem->stage_sub_prob[0].beta.size(); r++)
					problem->stage_sub_prob[0].beta[r].resize(problem->master_prob.num_var);
			}
			problem->solver.set_rhs_var(problem->stage_sub_prob[0]
				, *x, problem->stage_sub_prob[0].rhs, Cx);

			problem->solver.Solve_Prob(problem->stage_sub_prob[0], false);

			//fill out r_Cx
			int rv = 0;
			vec_d  r_Cx;
			vec_d  Cxx;
			for (int r = 0; r < problem->stage_sub_prob[0].rhs.size(); r++)
			{
				if (problem->stage_sub_prob[0].rhs[r].isRandom == true)
				{
					r_Cx.push_back(empt_scen.samp->at(rv) - problem->stage_sub_prob[0].Cx[r]);
					rv++;
				}
				Cxx.push_back(problem->stage_sub_prob[0].Cx[r]);
			}

			//sual of the stochastic constraints
			pi = problem->stage_sub_prob[0].dual_R;
			empt_clust.r_Cx = r_Cx;
			empt_clust.Cx = Cxx;
			Cxx.clear();
			r_Cx.clear();

			pi = problem->stage_sub_prob[0].dual_R;
			//scenario estimate
			empt_scen.est = problem->stage_sub_prob[0].zstar;
			//Cut Coeffs
			empt_scen.alpha = inner_product(problem->stage_sub_prob[0].dual_tot.begin(),
				problem->stage_sub_prob[0].dual_tot.end(), problem->stage_sub_prob[0].r_w.begin(), 0.0);

			int pos;
			empt_clust.new_pi = !in2Vector(*dual_pool, pi, pos);

			if (empt_clust.new_pi == true || iter == 0)
			{

				dual_pool->push_back(pi);
				dual_pool_num->push_back(1);
				empt_clust.pos = dual_pool_num->size() - 1;
				dual_num++;
				samples_clust->push_back(cluster_samples->size());

				//Sensitivity analysis ranges for the RHS
				/**empt_clust.SA_vals = problem->solver.getRHSSA_rand(problem->stage_sub_prob[0],
					*empt_clust.lower, *empt_clust.upper);
				*empt_clust.SA_z1 = sampl_ing->inRHSSA(problem->SPprobINFO.RVs, *empt_clust.SA_vals,
					*empt_clust.lower, *empt_clust.upper);
				empt_clust.z1 = empt_clust.SA_z1;*/
				create_z1(*empt_clust.z1, *empt_scen.z1);
				empt_clust.clust_size = 1.0;
				empt_clust.SSW = 0;

				//Cut Coeffs
				empt_scen.beta->resize(problem->master_prob.num_var);

				for (int r = 0; r < problem->stage_sub_prob[0].beta.size(); r++)
				{
					for (int v = 0; v < problem->master_prob.num_var; v++)
					{
						empt_scen.beta->at(v) += problem->stage_sub_prob[0].beta[r][v] * problem->stage_sub_prob[0].dual_tot[r];
					}
				}
				empt_clust.scens->push_back(empt_scen);

				//Estimation
				empt_clust.estim = empt_scen.est;
				empt_clust.alpha_estim = empt_scen.alpha;
				empt_clust.beta_estim = empt_scen.beta;
				empt_clust.estim_R = inner_product(pi.begin(), pi.end(), problem->stage_sub_prob[0].rho.begin(), 0.0);
				empt_clust.estim_D = empt_clust.estim - empt_clust.estim_R;
				//empt_clust.estim_lower = empt_clust.estim_D + inner_product(pi.begin(), pi.end(), empt_clust.lower->begin(), 0.0);
				//empt_clust.estim_upper = empt_clust.estim_D + inner_product(pi.begin(), pi.end(), empt_clust.upper->begin(), 0.0);
				//empt_clust.estim_extrm.push_back(empt_clust.estim_lower);
				//empt_clust.estim_extrm.push_back(empt_clust.estim_upper);
				//empt_clust.samples_extrm->push_back(*empt_clust.lower);
				//empt_clust.samples_extrm->push_back(*empt_clust.upper);
				empt_clust.SSW = 0;
				empt_clust.SSW = square_vec_d(empt_clust.estim_extrm);
				empt_clust.clust_size = 1;
				empt_clust.size = 1;

				if (print_SSW_info_seq == 1)
				{
					std::cout << "--------------------------------" << std::endl;
					std::cout << "cluster:  " << cluster_samples->size() + 1 << std::endl;
					print_vectord(empt_clust.samples->at(0));
					std::cout << "estim:  " << empt_clust.estim << std::endl;
					std::cout << "size:  " << empt_clust.clust_size << std::endl;
					std::cout << std::endl;
					std::cout << "dual: " << std::endl;
					print_vectord(pi);
					std::cout << std::endl << std::endl;
					print_vectord(*empt_clust.lower);
					std::cout << std::endl;
					print_vectord(*empt_clust.upper);
					std::cout << std::endl << std::endl;
					std::cout << "lower : " << empt_clust.estim_lower << "  -  upper : " <<
						empt_clust.estim_upper << std::endl;
					std::cout << "SSW:   " << empt_clust.SSW << std::endl;
					std::cout << "new pi:   " << std::boolalpha << empt_clust.new_pi << std::endl;
					if (empt_clust.new_pi == false)
					{
						//cout << dual_pool.size() << " " << empt_clust.pos << endl;
						print_vectord(dual_pool->at(empt_clust.pos));
						std::cout << std::endl;
						std::cout << dual_pool_num->at(empt_clust.pos) << std::endl;
					}
					std::cout << "--------------------------------" << std::endl;
				}

				*empt_clust.pi = pi;
				*empt_clust.pi_tot = problem->stage_sub_prob[0].dual_tot;

				//Add to the pool of clusters
				cluster_samples->push_back(empt_clust);
				scen_num++;

			}
			else
			{
				if (false) std::cout << scen_num << "  " << empt_scen.est << std::endl;
				if (cluster_samples->at(pos).active == true || iter < 100)
				{
					cluster_samples->at(pos).clust_size++;
					cluster_samples->at(pos).size++;
					update_z1(*cluster_samples->at(pos).z1, *empt_scen.z1);
					dual_pool_num->at(pos)++;
					empt_clust.pos = pos;
					add_scen2(problem, *empt_scen.samp, pos);
					samples_clust->push_back(pos);
					scen_num++;
				}
			}

			

		}
		else
		{
			int cl = samples_clust->at(samp_pos);
			cluster_samples->at(cl).clust_size++;
		}

		pi.clear();

		iter++;

		mean_MSW = 0;
		for (int cl = 0; cl < cluster_samples->size(); cl++)
			mean_MSW += cluster_samples->at(cl).SSW/(double) cluster_samples->size();
		if(false) std::cout << iter << "  "  << mean_MSW << std::endl;
		MSW_it.push_back(mean_MSW);
		if (iter > samp_num/(double)2) if (abs(MSW_it[iter - samp_num / (double)2] - mean_MSW)/ mean_MSW < 0.0001) break;
	}

	//system("pause");

	update_estimate_seq();

	Ll.clear();
	Ll.shrink_to_fit();
	Ul.clear();
	Ul.shrink_to_fit();
	pi.clear();
	pi.shrink_to_fit();
	Cx.clear();
	Cx.shrink_to_fit();

	if (print_SSW_info_seq == 1)
	{
		std::cout << estimate << std::endl;
		std::cout << tot_size << std::endl;
		int ss = getchar();
	}

	bool sampprint = false;
	if (sampprint == true)
	{
		for (int cl = 0; cl < cluster_samples->size(); cl++)
		{
			printf("before cl %d  -  size: %d  -  estim:  &%0.02f \n",
				cl, cluster_samples->at(cl).size, cluster_samples->at(cl).estim);
			print_vectord(*cluster_samples->at(cl).pi);
			for (int i = 0; i < cluster_samples->at(cl).scens->size(); i++)
				print_vectord(*cluster_samples->at(cl).scens->at(i).samp);
			print_vectord(*cluster_samples->at(cl).z1);
			printf("\n\n");
		}
	}

}
//**************************************************************************

//Get the size of the cluster by its SA range
double Point_Estim::getSize(Scen_cluster& clust)
{
	double size = 1;

	for (int i =0; i < clust.SA_vals->size(); i++)
	{
		size *= abs(clust.upper->at(i) - clust.lower->at(i));
	}

	return size;
}
//**************************************************************************

//**************************************************************************
// create neighbor N-one samples
void Point_Estim::create_neighbor_None(ProbPrep* problem, Sampling* sampl_ing, int clust, int num)
{

	vec2_d zero1 = *cluster_samples->at(clust).SA_z1;
	bool printing = false;


	for (int r = 0; r < problem->SPprobINFO.RVs.size(); r++)
	{
		vec_d new_r_samp1 = sampl_ing->real_in_range(problem->SPprobINFO.RVs[r], cluster_samples->at(clust).SA_vals->at(r));
		double min_est = new_r_samp1[0];
		double max_est = new_r_samp1[new_r_samp1.size() - 1];
		vec_d new_r_samp;
		new_r_samp.push_back(min_est);
		new_r_samp.push_back(max_est);

		for (int i = 0; i < new_r_samp.size(); i++)
		{

			vec_d samp = cluster_samples->at(clust).samples->at(0);
			samp[r] = new_r_samp[i];
			

			//The sample is new
			if (inVector2d(*cluster_samples->at(clust).samples, samp) == false)
			{
				Scen_struct  empt_scen;
				
				cluster_samples->at(clust).samples->push_back(samp);
				empt_scen.samp = &samp;
				empt_scen.active = true;
				samples->push_back(samp);
				samples_clust->push_back(clust);
				cluster_samples->at(clust).size++;

				//scenario estimate
				empt_scen.est = cluster_samples->at(clust).scens->at(0).est - 
					cluster_samples->at(clust).scens->at(0).samp->at(r) * cluster_samples->at(clust).pi->at(r) +
					samp[r] * cluster_samples->at(clust).pi->at(r);
				empt_scen.alpha = cluster_samples->at(clust).scens->at(0).alpha -
					              cluster_samples->at(clust).pi->at(r) * cluster_samples->at(clust).r_Cx[r] + 
					              (samp[r] - cluster_samples->at(clust).Cx[r]) * cluster_samples->at(clust).r_Cx[r];
				cluster_samples->at(clust).scens->push_back(empt_scen);

				if (printing == true)
				{
					printf("cluster:  %d   \n", clust);
					printf("current sample   \n");
					print_vectord(cluster_samples->at(clust).samples->at(0));
					printf("cluster size: %d   \n", cluster_samples->at(clust).size);
					std::cout << "new pi:  " << cluster_samples->at(clust).new_pi << std::endl;
					printf("pi val:   \n");
					print_vectord(*cluster_samples->at(clust).pi);
					std::cout << std::endl << std::endl;
					for (int rv = 0; rv < samp.size(); rv++)
					{
						printf("lower : %0.02f - new: %0.02f - main: %0.02f - upper: %0.02f  \n",
							cluster_samples->at(clust).SA_vals->at(rv)[0],
							samp[rv],
							cluster_samples->at(clust).samples->at(0).at(rv),
							cluster_samples->at(clust).SA_vals->at(rv)[1]);
					}
					printf("estimates   \n");
					printf("clust: %d  -  scen estim: %0.04f   -   cluster estim: %0.04f  \n", clust,
						empt_scen.est, cluster_samples->at(clust).estim);
					printf("estimate   \n");
					double est = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
						cluster_samples->at(clust).estim +
						(1 / (double)cluster_samples->at(clust).size) * empt_scen.est;
					std::cout << est << std::endl;
					printf("check its dual:   \n");
					print_vectord(problem->stage_sub_prob[0].dual_R);
					printf("---------------------------------------------------------------------- \n\n");
					int ss = getchar();
				}
				if (false)
				{
					printf("i(cluster): %d  -  r(rv): %d - w_i(r) - bar(w_i)(r): %0.02f  \n",
						clust, r, samp[r] - cluster_samples->at(clust).samples->at(0).at(r));
					int ss = getchar();
				}

				cluster_samples->at(clust).size_of_samps->push_back(1);

				scen_estim->push_back(empt_scen.est);
				cluster_samples->at(clust).scen_idx.push_back(scen_num);

				//Cut Coeffs
				cluster_samples->at(clust).alpha_estim =
					((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
					cluster_samples->at(clust).alpha_estim +
					(1 / (double)cluster_samples->at(clust).size) * empt_scen.alpha;

				//Estimation
				cluster_samples->at(clust).estim = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
					cluster_samples->at(clust).estim +
					(1 / (double)cluster_samples->at(clust).size) * empt_scen.est;

				scen_idx->at(clust).push_back(scen_num);
				scen_num++;
				
			}

			samp.clear();
			vec_d().swap(samp);
		}

		new_r_samp.clear();
		vec_d().swap(new_r_samp);

	}



}
//**************************************************************************

//**************************************************************************
// create neighbor N-one samples
void Point_Estim::create_neighbor_None_rmax(ProbPrep* problem, Sampling* sampl_ing, int clust, int num)
{

	vec2_d zero1 = *cluster_samples->at(clust).SA_z1;
	bool printing = false;
	int r = maxVector(*cluster_samples->at(clust).pi);

	vec_d new_r_samp = sampl_ing->real_in_range(problem->SPprobINFO.RVs[r], cluster_samples->at(clust).SA_vals->at(r));

	for (int i = 0; i < new_r_samp.size(); i++)
	{

		vec_d samp = cluster_samples->at(clust).samples->at(0);
		samp[r] = new_r_samp[i];


		//The sample is new
		if (inVector2d(*cluster_samples->at(clust).samples, samp) == false)
		{
			Scen_struct  empt_scen;

			cluster_samples->at(clust).samples->push_back(samp);
			empt_scen.samp = &samp;
			empt_scen.active = true;
			samples->push_back(samp);
			samples_clust->push_back(clust);
			cluster_samples->at(clust).size++;

			//scenario estimate
			empt_scen.est = cluster_samples->at(clust).scens->at(0).est -
				cluster_samples->at(clust).scens->at(0).samp->at(r) * cluster_samples->at(clust).pi->at(r) +
				samp[r] * cluster_samples->at(clust).pi->at(r);
			empt_scen.alpha = cluster_samples->at(clust).scens->at(0).alpha -
				cluster_samples->at(clust).pi->at(r) * cluster_samples->at(clust).r_Cx[r] +
				(samp[r] - cluster_samples->at(clust).Cx[r]) * cluster_samples->at(clust).r_Cx[r];
			cluster_samples->at(clust).scens->push_back(empt_scen);

			if (printing == true)
			{
				printf("cluster:  %d   \n", clust);
				printf("current sample   \n");
				print_vectord(cluster_samples->at(clust).samples->at(0));
				printf("cluster size: %d   \n", cluster_samples->at(clust).size);
				std::cout << "new pi:  " << cluster_samples->at(clust).new_pi << std::endl;
				printf("pi val:   \n");
				print_vectord(*cluster_samples->at(clust).pi);
				std::cout << std::endl << std::endl;
				for (int rv = 0; rv < samp.size(); rv++)
				{
					printf("lower : %0.02f - new: %0.02f - main: %0.02f - upper: %0.02f  \n",
						cluster_samples->at(clust).SA_vals->at(rv)[0],
						samp[rv],
						cluster_samples->at(clust).samples->at(0).at(rv),
						cluster_samples->at(clust).SA_vals->at(rv)[1]);
				}
				printf("estimates   \n");
				printf("clust: %d  -  scen estim: %0.04f   -   cluster estim: %0.04f  \n", clust,
					empt_scen.est, cluster_samples->at(clust).estim);
				printf("estimate   \n");
				double est = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
					cluster_samples->at(clust).estim +
					(1 / (double)cluster_samples->at(clust).size) * empt_scen.est;
				std::cout << est << std::endl;
				printf("check its dual:   \n");
				print_vectord(problem->stage_sub_prob[0].dual_R);
				printf("---------------------------------------------------------------------- \n\n");
				int ss = getchar();
			}
			if (true)
			{
				printf("i(cluster): %d  -  r(rv): %d - w_i(r) - bar(w_i)(r): %0.02f  \n",
					clust, r, samp[r] - cluster_samples->at(clust).samples->at(0).at(r));
				int ss = getchar();
			}

			cluster_samples->at(clust).size_of_samps->push_back(1);

			scen_estim->push_back(empt_scen.est);
			cluster_samples->at(clust).scen_idx.push_back(scen_num);

			//Cut Coeffs
			cluster_samples->at(clust).alpha_estim =
				((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
				cluster_samples->at(clust).alpha_estim +
				(1 / (double)cluster_samples->at(clust).size) * empt_scen.alpha;

			//Estimation
			cluster_samples->at(clust).estim = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
				cluster_samples->at(clust).estim +
				(1 / (double)cluster_samples->at(clust).size) * empt_scen.est;

			scen_idx->at(clust).push_back(scen_num);
			scen_num++;

		}

		samp.clear();
		vec_d().swap(samp);
	}

	new_r_samp.clear();
	vec_d().swap(new_r_samp);



}
//**************************************************************************
//**************************************************************************
void Point_Estim::create_neighbor_all(ProbPrep* problem, int clust)
{
	vec2_d tot_samp;
	int size__ = cluster_samples->at(clust).samples->size();
	for (int sc = 0; sc < size__; sc++)
	{
		vec_d orig_samp = cluster_samples->at(clust).samples->at(sc);
		Req_neighb_samp(tot_samp, orig_samp, 0, problem, clust);
		orig_samp.clear();
	}
}
//**************************************************************************
//**************************************************************************
//Generatig all neighbor samples
void Point_Estim::Req_neighb_samp(vec2_d& tot_samp, vec_d& orig_samp, int start, ProbPrep* problem, int clust)
{
	int end = problem->SPprobINFO.RVs.size();
	vec_d pi = *cluster_samples->at(clust).pi;

	if (end > start)
	{
		for (int rv = start; rv < end; rv++)
		{
			vec_d up_samp = neighb_samp(orig_samp, problem, rv, 1);
			if (inVector2d(*samples, up_samp) == false)
			{
				//tot_samp.push_back(up_samp);
				add_scen(problem, up_samp, clust);
			}	
			if (rv < end - 1)
			Req_neighb_samp(tot_samp, up_samp, rv + 1, problem,clust);

						
			vec_d dwn_samp = neighb_samp(orig_samp, problem, rv, 0);
			if (inVector2d(*samples, dwn_samp) == false)
			{
				//tot_samp.push_back(dwn_samp);
				add_scen(problem, dwn_samp, clust);
			}
			if (rv < end - 1)
			Req_neighb_samp(tot_samp, dwn_samp, rv + 1, problem,clust);


			if (inVector2d(*samples, orig_samp) == false)
			{
				//tot_samp.push_back(orig_samp);
				add_scen(problem, orig_samp, clust);
			}
			if (rv < end - 1)
			Req_neighb_samp(tot_samp, orig_samp, rv + 1, problem,clust);
		}
	}

}
//**************************************************************************
//**************************************************************************
//Generate neighbors of the given scenario
vec_d Point_Estim::neighb_samp(vec_d& orig_samp, ProbPrep* problem, int rv, int updwn)
{
	if (updwn == 0) //lower neighbor
	{
		RV_info rv_strct = problem->SPprobINFO.RVs[rv];
		std::map<double, double>::iterator it;
		vec_d output = orig_samp;
		
		double val = orig_samp[rv];

		for (it = rv_strct.CDF.begin(); it != rv_strct.CDF.end(); ++it)
		{
			if (val <= it->first)
			{
				break;
			}
			else
			{
				output[rv] = it->first;
			}
		}

		return output;

	}
	else if (updwn == 1) //upper neighborhood
	{
		RV_info rv_strct = problem->SPprobINFO.RVs[rv];
		std::map<double, double>::iterator it;
		vec_d output = orig_samp;

		double val = orig_samp[rv];

		for (it = rv_strct.CDF.begin(); it != rv_strct.CDF.end(); ++it)
		{
			if (val < it->first)
			{
				output[rv] = it->first;
				break;
			}
		}

		return output;

	}
	else
	{
		printf("ERROR: Point EStim: Neighborhood Generation \n");
		system("pause");
	}
}
//**************************************************************************

//**************************************************************************
//Add scenario to the cluster struct
void Point_Estim::add_scen(ProbPrep* problem, vec_d& samp, int clust)
{
	Scen_struct  empt_scen;
	empt_scen.samp = new vec_d;
	empt_scen.z1 = new vec_d;
	empt_scen.beta = new vec_d;
	*empt_scen.samp = samp;
	empt_scen.beta = cluster_samples->at(clust).beta_estim;

	//Put the sample in rhs struct
	int rv = 0;
	for (int r = 0; r < problem->stage_sub_prob[0].rhs.size(); r++)
	{
		if (problem->stage_sub_prob[0].rhs[r].isRandom == true)
		{
			problem->stage_sub_prob[0].rhs[r].value = (IloNum)empt_scen.samp->at(rv);
			rv++;
		}
	}

	problem->solver.set_rhs_var(problem->stage_sub_prob[0]
		, *x, problem->stage_sub_prob[0].rhs, problem->stage_sub_prob[0].Cx);
	problem->solver.Solve_Prob(problem->stage_sub_prob[0], false);
	empt_scen.est = problem->stage_sub_prob[0].zstar;
	
	if (false)
	{
		print_vectord(*cluster_samples->at(clust).pi);
		print_vectord(problem->stage_sub_prob[0].dual_R);
		print_vectord(*empt_scen.samp);
		std::cout << problem->stage_sub_prob[0].model << std::endl << std::endl;
		int ss = getchar();
	}


	if (problem->stage_sub_prob[0].dual_R == *cluster_samples->at(clust).pi)
	{
		samples->push_back(samp);
		empt_scen.alpha = inner_product(problem->stage_sub_prob[0].dual_tot.begin(),
			problem->stage_sub_prob[0].dual_tot.end(), problem->stage_sub_prob[0].r_w.begin(), 0.0);
		cluster_samples->at(clust).size++;
		cluster_samples->at(clust).scens->push_back(empt_scen);
		cluster_samples->at(clust).samples->push_back(samp);
		cluster_samples->at(clust).estim = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
			cluster_samples->at(clust).estim + (1 / (double)cluster_samples->at(clust).size) * empt_scen.est;
		cluster_samples->at(clust).alpha_estim = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
			cluster_samples->at(clust).alpha_estim + (1 / (double)cluster_samples->at(clust).size) * empt_scen.alpha;
		scen_num++;
	}
}
//**************************************************************************
//**************************************************************************
//Add scenario to the cluster struct for existing pi(pi is known)
void Point_Estim::add_scen2(ProbPrep* problem, vec_d& samp, int clust)
{
	Scen_struct  empt_scen;
	empt_scen.samp = new vec_d;
	empt_scen.beta = new vec_d;
	*empt_scen.samp = samp;
	empt_scen.beta = cluster_samples->at(clust).beta_estim;

	//Put the sample in rhs struct
	int rv = 0;
	vec_d r_Cx;
	for (int r = 0; r < problem->stage_sub_prob[0].rhs.size(); r++)
	{
		if (problem->stage_sub_prob[0].rhs[r].isRandom == true)
		{
			r_Cx.push_back(samp[rv] - problem->stage_sub_prob[0].Cx[r]);
			rv++;
		}
	}

	
	empt_scen.est = cluster_samples->at(clust).estim_D + inner_product(cluster_samples->at(clust).pi->begin(),
		cluster_samples->at(clust).pi->end(), r_Cx.begin(), 0.0);;

	if (false)
	{
		print_vectord(*cluster_samples->at(clust).pi);
		print_vectord(problem->stage_sub_prob[0].dual_R);
		print_vectord(*empt_scen.samp);
		std::cout << problem->stage_sub_prob[0].model << std::endl << std::endl;
		int ss = getchar();
	}


	empt_scen.alpha = inner_product(problem->stage_sub_prob[0].dual_tot.begin(),
		problem->stage_sub_prob[0].dual_tot.end(), problem->stage_sub_prob[0].r_w.begin(), 0.0);
	cluster_samples->at(clust).scens->push_back(empt_scen);
	cluster_samples->at(clust).samples->push_back(samp);
	cluster_samples->at(clust).estim = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
		cluster_samples->at(clust).estim + (1 / (double)cluster_samples->at(clust).size) * empt_scen.est;
	cluster_samples->at(clust).alpha_estim = ((cluster_samples->at(clust).size - 1) / (double)cluster_samples->at(clust).size) *
		cluster_samples->at(clust).alpha_estim + (1 / (double)cluster_samples->at(clust).size) * empt_scen.alpha;
	double SSW = 0;
	for (int j = 0; j < cluster_samples->at(clust).scens->size(); j++)
		SSW += (cluster_samples->at(clust).scens->at(j).est - cluster_samples->at(clust).estim)*
		(cluster_samples->at(clust).scens->at(j).est - cluster_samples->at(clust).estim);
	cluster_samples->at(clust).SSW_it.push_back(SSW);
	if (std_vec_d(cluster_samples->at(clust).SSW_it)/ cluster_samples->at(clust).estim  > 0.01)
	{
		if (cluster_samples->at(clust).SSW_it.size() > 100)
		{
			if (false)
			{
				printf("cl: %d - size: %d -  %0.02f  \n", clust, cluster_samples->at(clust).SSW_it.size(),
					abs(cluster_samples->at(clust).SSW_it[cluster_samples->at(clust).SSW_it.size() - 100] - SSW) / SSW);
			}
			if (abs(cluster_samples->at(clust).SSW_it[cluster_samples->at(clust).SSW_it.size() - 100] - SSW)/SSW < 0.01)
			{
				cluster_samples->at(clust).active = false;
			}
		}
	}

}
//**************************************************************************

//**************************************************************************
//Create the z1 of the cluster
void Point_Estim::create_z1(vec2_d& z1_main, vec_d& z1)
{
	z1_main.resize(z1.size());
	for (int r = 0; r < z1.size(); r++)
	{
		z1_main[r].push_back(std::max(0.0,z1[r] - 0.1));
		z1_main[r].push_back(std::min(1.0,z1[r] + 0.1));
		//z1_main[r].push_back(0.0);
		//z1_main[r].push_back(1.0);
	}
}
//**************************************************************************

//**************************************************************************
//update the z1 of the cluster
void Point_Estim::update_z1(vec2_d& z1_main, vec_d& z1)
{
	for (int r = 0; r < z1.size(); r++)
	{
		z1_main[r][0] = std::max(0.0,std::min(z1_main[r][0], z1[0]) - 0.1);
		z1_main[r][1] = std::min(1.0,std::max(z1_main[r][1], z1[1]) + 0.1);
	}
}
//**************************************************************************
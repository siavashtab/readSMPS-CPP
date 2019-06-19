#include "LShaped_Algo.h"



LShaped_Algo::LShaped_Algo()
{
	std::cout << "--------------------------------------------------------------------------------" << std::endl;

	std::cout << "                               LShaped Algorithm                                " << std::endl;

	std::cout << "--------------------------------------------------------------------------------" << std::endl;
	std::cout << "--------------------------------------------------------------------------------" << std::endl;
}


LShaped_Algo::~LShaped_Algo()
{
	//cout << "--------------------------------------------------------------------------------" << endl;
	std::cout << "------------------------------ LShaped Class is Closed -------------------------" << std::endl;
	Beta.clear();
	Beta.shrink_to_fit();
	Alpha_s.clear();
	Alpha_s.shrink_to_fit();
	Beta_s.clear();
	Beta_s.shrink_to_fit();
	pi.clear();
	pi.shrink_to_fit();
	pi_collect.clear();
	pi_collect.shrink_to_fit();
	LShapedINFO.GAP_it.clear();
	LShapedINFO.GAP_it.shrink_to_fit();
	LShapedINFO.LB_iter.clear();
	LShapedINFO.LB_iter.shrink_to_fit();
	LShapedINFO.UB_now.clear();
	LShapedINFO.UB_now.shrink_to_fit();
	LShapedINFO.UB_w_iter.clear();
	LShapedINFO.UB_w_iter.shrink_to_fit();
	LShapedINFO.SSW_vec.clear();
	LShapedINFO.SSW_vec.shrink_to_fit();
	LShapedINFO.xFirst.clear();
	LShapedINFO.xFirst.shrink_to_fit();
	LShapedINFO.x_star.clear();
	LShapedINFO.xFirst.shrink_to_fit();
	cut_type_it.clear();
	cut_type_it.shrink_to_fit();
	LShaped_cut_tot.clear();
	LShaped_cut_tot.shrink_to_fit();
	LShaped_cut.clear();
	LShaped_cut.shrink_to_fit();
	LShapedINFO.GAP = 100;
	LShapedINFO.LB = -INFINITY;
	LShapedINFO.UB = INFINITY;
	samp_clust.clear();
	samp_clust.shrink_to_fit();
	scen_pi_norm.clear();
	scen_pi_norm.shrink_to_fit();
}

void LShaped_Algo::Initialize(ProbPrep* problem)
{
	cut_num = 0;
	nScenario = problem->SPprobINFO.initial_scen_num;
	LShapedINFO.Algo_isStable = false;
	LShapedINFO.GAP = 100;
	LShapedINFO.LB = -INFINITY;
	LShapedINFO.UB = INFINITY;
	LShapedINFO.UB_now.clear();
	LShapedINFO.UB_w_iter.clear();
	LShapedINFO.max_iter = LShaped_maxiter;
	LShapedINFO.min_iter = LShaped_miniter;
	LShapedINFO.sigma = LShaped_reg_sigma;
	for (int v = 0; v < problem->master_prob.num_var; v++)
	{
		LShapedINFO.xFirst.push_back(problem->mean_prob.sol[v]);
	}	
	LShapedINFO.x_star = LShapedINFO.xFirst;
	LShapedINFO.opt_pi.resize(nScenario);
	iter = 0;
	Optimality = false;
	feasibility = false;
}

void LShaped_Algo::optimize(ProbPrep* problem, Sampling* sampl_ing, std::vector<Partitions>& parts, int siter)
{
	Initialize(problem);

	//initialize the point estimation class
	p_estim = new Point_Estim();
	p_estim->initialize(nScenario, siter);
	cumul_tot_scen = 0; 

	std::string file_name_LShaped = problem->SPprobINFO.dirr_algo +
		                       problem->SPprobINFO.output_algo +
		                       std::to_string(siter) +
		                       problem->SPprobINFO.algo_type;

	if(print_Algo_results == 1) std::ofstream file_LShaped(file_name_LShaped);
	problem->master_prob.zstar_wout_surro = abs(problem->mean_prob.zstar);
	
	//initialize running time computation
	opt_time = 0;
	clock_t topt1,topt2;
	topt1 = clock();

	//Main Algorithm Loop
	while (stopping_check())
	{
		if (iter % 100 == 0) std::cout << std::endl;
		printf("*");
		
		
		if(iter > 0) solve_master(problem);

		if (LShapedINFO.LB > LShapedINFO.UB)
		{
			LShapedINFO.LB = LShapedINFO.UB;
			break;
		}
		
		solve_sub(problem, sampl_ing, parts);
		
		if (display_LShaped_LB_UB == 1)
		{
			printf(" LB: %0.02f   -   UB: %0.02f  , ", LShapedINFO.LB, LShapedINFO.UB);
		}

		LShapedINFO.GAP = (LShapedINFO.UB - LShapedINFO.LB) * 100 / abs(LShapedINFO.UB);



		iter++;
	}
	
	//capture the optimizaion time
	topt2 = clock();
	opt_time = ((double)topt2 - (double)topt1) / CLOCKS_PER_SEC;

	p_estim->end();

	//displaying the point estimation
	if (disp_point_estim == 1)
	{
		Point_Estimation_Analysis(problem, sampl_ing);
	}


	if (print_Master_At_end == 1)
	{
		std::cout << problem->master_prob.model << std::endl << std::endl;
		printf("LShaped Cuts:    \n\n");
		for (int i = 0; i < problem->master_prob.LShaped_opt.getSize(); i++)
		{
			printf("Cut Number :  %d  \n", i);
			std::cout << problem->master_prob.LShaped_opt[i] << std::endl;
			printf("press a key to continue    \n\n");
			std::string nothing;
			std::cin >> nothing;
		}
		system("pause");
	}



	//Close generated cuts
	if (Optimality == true)
	{
		printf("c");
		problem->solver.Clear_Prob(problem->master_prob, "optimality");
		printf("c");
	}
	
	if (feasibility == true)
	{
		printf("c");
		problem->solver.Clear_Prob(problem->master_prob, "feasibility");
		printf("c");
	}

	p_estim->end();
	delete p_estim;
	//system("pause");

}

//This subrouting is just for point estimation analysis with different sampling techniques
void LShaped_Algo::Point_Estimation_Analysis(ProbPrep* problem, Sampling* sampl_ing)
{
	printf("\n");
	
	//initialize running time computation
	double optt_time = 0;
	clock_t topt1, topt2;

	//Optimal Point Estimation
	p_estim->x = &LShapedINFO.xFirst;


	printf("        ,    ,    Monte     ,          ,");
	if (scen_or_clust == 0) if(scen_type == 1)  printf("        ,    ,    LHS       ,          ,");
	printf("        ,    ,    ACS       ,          ,\n");
	printf("scen  ,   est   ,  std   ,  time  , duals ,  ");
	if (scen_or_clust == 0) if (scen_type == 1)  printf("   est   ,  std   ,  time  , duals , ");
	printf("   est   ,  std   ,  time  , size  ,  duals \n");
	printf("-------------------------------------");
	printf("-------------------------------------");
	printf("-------------------------------------\n");

	for (int s = 1; s < 30; s++)
	{
		int scen;
		if(s <= 10) scen = s * 10;
		if(s > 10) scen = (s-10 + 1) * 100;
		if (s > 20) scen = (s - 20 + 1) * 1000;

		vec_d estim1;
		vec_d estim2;
		vec_d estim3;
		vec_d tim1;
		vec_d tim2;
		vec_d tim3;
		vec_d scen_size;
		vec_d ratio;
		vec_d dual_1;
		vec_d dual_2;
		vec_d dual_3;

		for (int rep = 0; rep < 30; rep++)
		{
			topt1 = clock();
			p_estim->initialize(scen, rep);
			p_estim->monte(problem, sampl_ing);
			topt2 = clock();
			optt_time = ((double)topt2 - (double)topt1) / CLOCKS_PER_SEC;
			estim1.push_back(p_estim->estimate);
			tim1.push_back(optt_time);
			if(false) std::cout << p_estim->estimate << " " << p_estim->scen_num << std::endl;
			dual_1.push_back((double)p_estim->dual_num);
			p_estim->end();
			
			if (scen_or_clust == 0)
			{
				if (scen_type == 1)
				{
					topt1 = clock();
					p_estim->lhs(problem, sampl_ing);
					topt2 = clock();
					optt_time = ((double)topt2 - (double)topt1) / CLOCKS_PER_SEC;
					estim2.push_back(p_estim->estimate);
					tim2.push_back(optt_time);
					if (false) std::cout << p_estim->estimate << " " << p_estim->scen_num << std::endl;
					dual_2.push_back((double)p_estim->dual_num);
					p_estim->end();
				}
			}


			topt1 = clock();
			p_estim->acs(problem, sampl_ing);
			topt2 = clock();
			optt_time = ((double)topt2 - (double)topt1) / CLOCKS_PER_SEC;
			estim3.push_back(p_estim->estimate);
			tim3.push_back(optt_time);
			scen_size.push_back((double)p_estim->scen_num);
			if (false)std::cout << p_estim->estimate << " " << p_estim->scen_num << std::endl;
			double r2 = p_estim->R_a;
			ratio.push_back(r2);
			dual_3.push_back((double)p_estim->dual_num);
			if (false)
			{
				if (scen == 10)
				{
					
					for (int i = 0; i < p_estim->cluster_samples->size(); i++)
					{
						printf("cluster  %d  --  estim:  %0.02f  -- estim D: %0.02f \n", i, p_estim->cluster_samples->at(i).estim,
							p_estim->cluster_samples->at(i).estim_D);
						printf("Cx: \n");
						print_vectord(p_estim->cluster_samples->at(i).Cx);
						printf("dual: \n");
						print_vectord(*p_estim->cluster_samples->at(i).pi);
						printf("samples size:  %d  - samples: \n", p_estim->cluster_samples->at(i).samples->size());
						for (int j = 0; j < p_estim->cluster_samples->at(i).samples->size(); j++)
						{
							printf("%d  - est: %0.02f  -  ", j , p_estim->cluster_samples->at(i).scens->at(j).est);
							print_vectord(*p_estim->cluster_samples->at(i).scens->at(j).samp);
						}
						printf("\n\n");
					}
					int ss = getchar();
				}
			}
			p_estim->end();
			
		}
		
		if (scen_or_clust == 0 && scen_type == 1)
		{
			if (std_vec_d(estim1) < std_vec_d(estim3) && std_vec_d(estim1) < std_vec_d(estim2))
			{
				printf("%d   ,   %0.2f  , * %0.4f , %0.2f  , %0.2f , ", scen,
					mean_vec_d(estim1), std_vec_d(estim1), mean_vec_d(tim1), mean_vec_d(dual_1));
			}
			else
			{
				printf("%d   ,   %0.2f  ,  %0.4f , %0.2f  , %0.2f  , ", scen,
					mean_vec_d(estim1), std_vec_d(estim1), mean_vec_d(tim1), mean_vec_d(dual_1));
			}

			if (std_vec_d(estim2) < std_vec_d(estim1) && std_vec_d(estim2) < std_vec_d(estim3))
			{
				printf(" %0.2f  , * %0.4f , %0.2f  , %0.2f , ",
					mean_vec_d(estim2), std_vec_d(estim2), mean_vec_d(tim2), mean_vec_d(dual_2));
			}
			else
			{
				printf(" %0.2f  ,  %0.4f , %0.2f  , %0.2f , ",
					mean_vec_d(estim2), std_vec_d(estim2), mean_vec_d(tim2), mean_vec_d(dual_2));
			}
		}
		else
		{
			if (std_vec_d(estim1) < std_vec_d(estim3))
			{
				printf("%d   ,   %0.2f  , * %0.4f , %0.2f  , %0.2f , ", scen,
					mean_vec_d(estim1), std_vec_d(estim1), mean_vec_d(tim1), mean_vec_d(dual_1));
			}
			else
			{
				printf("%d   ,   %0.2f  ,  %0.4f , %0.2f  , %0.2f  , ", scen,
					mean_vec_d(estim1), std_vec_d(estim1), mean_vec_d(tim1), mean_vec_d(dual_1));
			}
		}

		if (scen_or_clust == 0 && scen_type == 1)
		{
			if (std_vec_d(estim3) < std_vec_d(estim1) && std_vec_d(estim3) < std_vec_d(estim2))
			{
				printf("%0.2f  , * %0.4f , %0.2f , %0.2f , %0.2f \n",
					mean_vec_d(estim3), std_vec_d(estim3), mean_vec_d(tim3),
					mean_vec_d(scen_size), mean_vec_d(dual_3));
			}
			else
			{
				printf("%0.2f  ,  %0.4f , %0.2f , %0.2f , %0.2f , \n",
					mean_vec_d(estim3), std_vec_d(estim3), mean_vec_d(tim3),
					mean_vec_d(scen_size), mean_vec_d(dual_3));
			}
		}
		else
		{
			if (std_vec_d(estim3) < std_vec_d(estim1))
			{
				printf("%0.2f  , * %0.4f , %0.2f , %0.2f , %0.2f \n",
					mean_vec_d(estim3), std_vec_d(estim3), mean_vec_d(tim3),
					mean_vec_d(scen_size), mean_vec_d(dual_3));
			}
			else
			{
				printf("%0.2f  ,  %0.4f , %0.2f , %0.2f , %0.2f , \n",
					mean_vec_d(estim3), std_vec_d(estim3), mean_vec_d(tim3),
					mean_vec_d(scen_size), mean_vec_d(dual_3));
			}
		}

		
	}

	system("pause");

}

bool LShaped_Algo::stopping_check()
{
	if (iter > LShapedINFO.min_iter)
	{
		if (LShapedINFO.GAP < LShaped_gap)
		{
			return false;
		}
		else if (iter > LShapedINFO.max_iter)
		{
			return false;
		}
		else
		{
			return true;
		}
	}
	else
	{
		return true;
	}
}

void LShaped_Algo::solve_sub(ProbPrep* problem, Sampling* sampl_ing, std::vector<Partitions>& parts)
{
	nScenario = problem->SPprobINFO.initial_scen_num;
	LShaped_cut_tot.push_back(LShaped_cut);
	LShaped_cut.clear();	
	LShapedINFO.UB_now.clear();
	double sum = 0;
	double ub = problem->master_prob.zstar_wout_surro;

	if (Sampling_Method == 0)
	{
		p_estim->x = &LShapedINFO.xFirst;
		try
		{
			p_estim->monte(problem, sampl_ing);
			if (disp_point_estim_iter == 1) Point_Estimation_Analysis(problem, sampl_ing);
			cumul_tot_scen += p_estim->scen_num;
		}
		catch (...)
		{
			printf("ERROR: Lshaped - solve_sub: point estim: monte!");
			system("pause");
		}		
		ub += p_estim->estimate;
		create_cuts(problem);
		put_opt_in_range(problem);
	}
	else if (Sampling_Method == 1)
	{
		p_estim->x = &LShapedINFO.xFirst;
		try
		{
			p_estim->lhs(problem, sampl_ing);
			if(disp_point_estim_iter == 1) Point_Estimation_Analysis(problem, sampl_ing);
			cumul_tot_scen += p_estim->scen_num;
		}
		catch (...)
		{
			printf("ERROR: Lshaped - solve_sub: point estim: lhs!");
			system("pause");
		}		
		ub += p_estim->estimate;
		create_cuts(problem);
		put_opt_in_range(problem);
	}
	else
	{
		p_estim->x = &LShapedINFO.xFirst;
		if (LShapedINFO.UB_w_iter.size() > 1 && LShapedINFO.UB_w_iter[LShapedINFO.UB_w_iter.size() - 1] ==
			LShapedINFO.UB_w_iter[LShapedINFO.UB_w_iter.size() - 2])
		{
			try
			{
				p_estim->lhs(problem, sampl_ing);
				if (disp_point_estim_iter == 1) Point_Estimation_Analysis(problem, sampl_ing);
				cumul_tot_scen += p_estim->scen_num;
			}
			catch (...)
			{
				printf("ERROR: Lshaped - solve_sub: point estim: lhs!");
				system("pause");
			}
		}
		else
		{
			try
			{
				p_estim->acs(problem, sampl_ing);
				if (disp_point_estim_iter == 1) Point_Estimation_Analysis(problem, sampl_ing);
				cumul_tot_scen += p_estim->scen_num;
			}
			catch (...)
			{
				printf("ERROR: Lshaped - solve_sub: point estim: acs!");
				system("pause");
			}
		}	
		ub += p_estim->estimate;
		create_cuts(problem);
		put_opt_in_range(problem);
	}

	if ((ub< 0.7*LShapedINFO.UB) && (Regular_LShaped == 1))
	{
		LShapedINFO.x_star = LShapedINFO.xFirst;
		problem->master_prob.update_xstar = true;
	}

	LShapedINFO.UB = ub;
	LShapedINFO.UB_w_iter.push_back(ub);

	if (show_master_iter == 1)
	{
		std::cout << problem->master_prob.model << std::endl;
		std::string nothing;
		std::cin >> nothing;
	}
	
	if (display_estimator_iter0 == 1 && iter == 0)
	{
		LShapedINFO.UB_iter = LShapedINFO.UB;
		LShapedINFO.UB_std_iter = std_vec_d(LShapedINFO.UB_now);
	}

	
}


#pragma region Opt LShaped cut preparation functions
void LShaped_Algo::create_cuts(ProbPrep* problem)
{

	Optimality = true;

	//Creating the current set of cuts
	for (int i = 0; i < p_estim->cluster_samples->size(); i++)
	{
		//empty struct
		LShaped_Cut empty_cut;
		empty_cut.weight = 1 / (double)(p_estim->cluster_samples->size());
		double weight = (1 / (double)(p_estim->cluster_samples->size()));
		empty_cut.Alpha_s = weight * p_estim->cluster_samples->at(i).alpha_estim;
		empty_cut.Beta_s = const_prod(*p_estim->cluster_samples->at(i).beta_estim, weight);

		empty_cut.cut_type = "optimality";
		empty_cut.iter = iter;
		empty_cut.cut_num = cut_num;
		empty_cut.multp = 1;
		cut_num++;

		//Creating the cut
		IloExpr expr(*problem->master_prob.env);
		expr = problem->master_prob.surro_vars_raw[i];
		for (int v = 0; v < problem->master_prob.num_var; v++) expr += 
			*LShapedINFO.xFirst[v].var * empty_cut.Beta_s[v];
		empty_cut.cut = expr >= empty_cut.Alpha_s;
		
		char varname[20];
		sprintf(varname, "\cut_%d_%d", (int)(iter), (int)(i));
		empty_cut.cut.setName(varname);
		empty_cut.Name = varname;
		expr.clear();

		//Put it in the pool
		LShaped_cut.push_back(empty_cut);
	}

	p_estim->end();

}
void LShaped_Algo::put_opt_in_range(ProbPrep* problem)
{

	problem->solver.open_rng(problem->master_prob.LShaped_opt_iter, iter);
	problem->solver.open_rng(problem->master_prob.LShaped_opt, iter);
	//cout << "model created by LShaped cuts!" << endl;

	for (int i = 0; i < LShaped_cut.size(); i++)
	{
		problem->master_prob.LShaped_opt_iter.add(LShaped_cut[i].cut);
	}

	problem->master_prob.LShaped_opt.add(problem->master_prob.LShaped_opt_iter);


	
}
#pragma endregion Opt LShaped cut preparation functions

void LShaped_Algo::solve_master(ProbPrep* problem)
{
		
	if (Regular_LShaped == 1)
	{
		problem->x_reg = LShapedINFO.x_star;
		problem->master_prob.xhat = LShapedINFO.x_star;
		problem->master_prob.isReg = true;
		if (problem->master_prob.update_xstar)
		{
			problem->master_prob.model.remove(problem->master_prob.obj_raw);
			problem->add_master_obj(problem->master_prob, problem->solver);
			problem->master_prob.model.add(problem->master_prob.obj_raw);
		}
	}
	
	problem->master_prob.model.add(problem->master_prob.LShaped_opt_iter);

	
	//problem.master_prob.cplex.exportModel("model.lp");
	problem->solver.Solve_Prob(problem->master_prob, false);
	

	LShapedINFO.xFirst = problem->master_prob.sol;
	vec_d sol;
	for (int i = 0; i < problem->master_prob.sol.size(); i++)
		sol.push_back(problem->master_prob.sol[i].value);

	LShapedINFO.LB = problem->master_prob.zstar;
	LShapedINFO.LB_iter.push_back(LShapedINFO.LB);

	if (display_Master_surrogates == 1)
	{
		std::cout << std::endl;
		std::cout << problem->master_prob.model << std::endl;
		std::cout << "values: " << std::endl;
		for (int i = 0; i < problem->master_prob.sol_surrogate.size(); i++)
		{
			printf(" %s  :  %0.2f   ,   ", problem->master_prob.sol_surrogate[i].col, problem->master_prob.sol_surrogate[i].value);
			if (i % 4 == 0) std::cout << std::endl;
		}
		int empty;
		std::cin >> empty;
	}
	
	problem->master_prob.LShaped_opt_iter.clear();


	//enhance the QP convergence
	if (iter > LShapedINFO.max_iter / (double)2)
	{
		problem->master_prob.sigma += 0.01;
	}

	

}



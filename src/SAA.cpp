
#include "SAA.h"
#include"seed_SAA.h"

SAA::SAA()
{

	printf("--------------------------------------------------------------------------------\n");
	printf("                          SAMPLE AVERAGE APPROXIMATION                          \n");
	printf("--------------------------------------------------------------------------------\n");
	printf("--------------------------------------------------------------------------------\n");
	
}

void SAA::reset_SAA()
{
	printf("--------------------------------------------------------------------------------\n");
	printf("----------------------------- SAA Class is Closed ------------------------------\n");
	rep_UB.clear();
	rep_UB.shrink_to_fit();
	rep_LB.clear();
	rep_LB.shrink_to_fit();
	rep_partitions.clear();
	rep_partitions.shrink_to_fit();
	rep_Scen.clear();
	rep_Scen.shrink_to_fit();
	parts.clear();
	parts.shrink_to_fit();
	opt_time = 0;
}

void SAA::Initialize(ProbPrep*  problem)
{
	N = problem->SPprobINFO.initial_scen_num;
	part_num = 0;
	opt_time = 0;

	try
	{
		optimize(problem);
	}
	catch (const std::exception&)
	{
		std::cerr << "SAA: Initialize: Optimize!" << std::endl;
		system("pause");
	}

}

void SAA::optimize(ProbPrep*  problem)
{

	clock_t t1, t2;
	t1 = clock();
	double seconds1 = 0;

	std::string file_name_saa = problem->SPprobINFO.dirr_saa +
		problem->SPprobINFO.output_saa +
		std::to_string(problem->SPprobINFO.initial_scen_num) +
		problem->SPprobINFO.algo_type;
	//ofstream file_saa(file_name_saa);


	//write_legend(file_saa, problem.SPprobINFO.file_name);

	//cout << "--------------------------------------------------------------------------------" << endl;

	
	for (int siter = 0; siter < SAA_rep_num; siter++)
	{

		LShaped_Algo*   LShaped;
		LShaped = new LShaped_Algo();
		Sampling* sampl_ing;
		sampl_ing = new Sampling();

		try
		{
			LShaped->optimize(problem, sampl_ing, parts, siter);
			part_num = LShaped->pi_collect.size();
		}
		catch (const std::exception&)
		{
			std::cerr << "SAA: SAMPLING: Algo Optimize!" << std::endl;
			system("pause");
		}

		if (display_estimator_iter0 == 1) UB_iter_firstIter.push_back(LShaped->LShapedINFO.UB_iter);

		opt_time += LShaped->opt_time;

		x = LShaped->LShapedINFO.xFirst;
		int ub_it = eval(problem, sampl_ing, siter);

		t2 = clock();
		double diff1((double)t2 - (double)t1);
		seconds1 = diff1 / CLOCKS_PER_SEC;

		rep_LB.push_back(LShaped->LShapedINFO.LB);
		rep_partitions.push_back(part_num);
		rep_Scen.push_back(LShaped->cumul_tot_scen);

		printf("\n");
		printf("--------------------------------------------------------------------------------\n");

		if (Regular_LShaped == 1) {
			printf("Algorithm: Regularized L-Shaped \n");
		}
		else if (Regular_LShaped == 0) {
			printf("Algorithm: L-Shaped \n");
		}
		if (Sampling_Method == 0) {
			printf("Sampling: Monte Carlo \n");
		}
		else if (Sampling_Method == 1) {
			printf("Sampling: Latin Hypercube Sampling \n");
		}

		printf("\n");
		printf("Rep        :%d \n", siter + 1);
		printf("LB         :%0.2f \n", rep_LB[siter]);
		printf("UB         :%0.2f -- # of Evaluated Scenarios: %d -- Eval STD: %0.2f \n", rep_UB[siter], ub_it, stdes);
		printf("Parts      :%d \n", part_num);
		printf("Scen       :%d \n", sampl_ing->samples.size());
		printf("Cut num    :%d \n", LShaped->cut_num);
		printf("opt time   :%0.4f \n", opt_time);
		printf("time       :%0.3f \n", seconds1);
		printf("--------------------------------------------------------------------------------\n");

		if(false)
		{
			//file_saa << "--------------------------------------------------------------------------------" << endl;
			//file_saa << "Rep  ," << siter + 1 << ", LB =  ," << rep_LB[siter] << ", UB =  ," << rep_UB[siter]
			//	<< ", STD. = ," << stdes;
			//file_saa << ", -- UB Eval Iterations = ,"<< ub_it ;
			//file_saa << ", Parts # =  ," << part_num;
			//file_saa << ", Scen # =  ," << sampl_ing.samples.size();
			//file_saa << ", opt t =  ," << opt_time;
			//file_saa << ", t =  ," << seconds1 << endl;
			//file_saa << "--------------------------------------------------------------------------------" << endl;		
		}


		//sampl_ing->~Sampling();
		//LShaped->~LShaped_Algo();
		delete sampl_ing;
		delete LShaped;

	}

	double sqrtM = std::sqrt(SAA_rep_num);
	double LBsum = std::accumulate(rep_LB.begin(), rep_LB.end(), 0.0);
	double LBmean = LBsum / (double)rep_LB.size();
	double LB_std = std_vec_d(rep_LB);
	printf("\n");
	printf("LB mean: %0.2f  -- LB Std.: %0.2f \n", LBmean, LB_std);
	printf("LB CI:  [ %0.2f , %0.2f] \n\n", LBmean - 1.96 * LB_std / sqrtM, LBmean + 1.96 * LB_std / sqrtM);

	double UBsum = std::accumulate(rep_UB.begin(), rep_UB.end(), 0.0);
	double UBmean = UBsum / (double)rep_UB.size();
	double UB_std = std_vec_d(rep_UB);
	printf("UB mean: %0.2f  -- UB Std.: %0.2f \n", UBmean, UB_std);
	printf("UB CI:  [ %0.2f , %0.2f] \n\n", UBmean - 1.96 * UB_std / sqrtM, UBmean + 1.96 * UB_std / sqrtM);

	printf("Pessimistic Gap: %0.2f \n", UBmean + 1.96 * UB_std / sqrtM - (LBmean - 1.96 * LB_std / sqrtM));

	int partssum = std::accumulate(rep_partitions.begin(), rep_partitions.end(), 0.0);
	double partmean = partssum / (double)rep_partitions.size();
	double part_std = 1.96 * std_vec_i(rep_partitions) / sqrtM;

	int scensum = std::accumulate(rep_Scen.begin(), rep_Scen.end(), 0.0);
	double scenmean = scensum / (double)rep_Scen.size();
	double scen_std = 1.96 * std_vec_i(rep_Scen) / sqrtM;

	AlgOutput.LB_mean = LBmean;
	AlgOutput.LB_std  = LB_std;
	AlgOutput.LB_CI_std = 1.96 * LB_std / (double)sqrtM;
	AlgOutput.UB_mean = UBmean;
	AlgOutput.UB_std = UB_std;
	AlgOutput.UB_CI_std = 1.96 * UB_std / (double)sqrtM;
	AlgOutput.Pessim_gap = UBmean + 1.96 * UB_std / (double)sqrtM - (LBmean - 1.96 * LB_std / (double)sqrtM);
	AlgOutput.mean_part_num = part_num;
	AlgOutput.opt_time = opt_time;
	AlgOutput.mean_part_num = partmean;
	AlgOutput.std_part_num = part_std;
	AlgOutput.mean_scen_num = scenmean;
	AlgOutput.std_scen_num = scen_std;

	if (false)
	{
		//file_saa << "LB mean = ," << LBmean << " +- " << 1.96 * LB_std / sqrtM << " , ";
		//file_saa << "LB CI: ," << "[ ," << LBmean - 1.96 * LB_std << " , " << LBmean + 1.96 * LB_std << ", ]" << endl << endl;
		//file_saa << "UB mean = ," << UBmean << " +- " << 1.96 * UB_std / sqrtM << " , ";
		//file_saa << "UB CI: ," << "[ ," << UBmean - 1.96 * UB_std << " , " << UBmean + 1.96 * UB_std << ", ]" << endl << endl;
		//file_saa << "Cluster mean = ," << mean_vec_d(rep_partitions) << " +- " << 1.96 * std_vec_d(rep_partitions) / sqrtM << " , ";
		//file_saa << "Scen mean = ," << mean_vec_d(rep_Scen) << " +- " << 1.96 * std_vec_d(rep_Scen) / sqrtM << " , ";
		//file_saa << "Mean Gap:  ," << UBmean - LBmean << " , ";
		//file_saa << "Pessimistic Gap:  ," << UBmean + 1.96 * UB_std / sqrtM - (LBmean - 1.96 * LB_std / sqrtM) << " , ";
		//file_saa << "OPT time:  ," << opt_time << " , ";
		//file_saa << "TOT time:  ," << seconds1 << endl << endl;
	}


	if (display_estimator_iter0 == 1)
	{
		printf("The variablility of the UBs of iter 0: %0.2f \n", std_vec_d(UB_iter_firstIter));
		//file_saa << "The variablility of the UBs of iter 0:  " << std_vec_d(UB_iter_firstIter) << endl << endl;
	}


	reset_SAA();
}


int SAA::eval(ProbPrep* problem, Sampling* sampl_ing, int siter)
{
	printf(">>");

	double UB = 0;
	stdes = INFINITY;
	double UBmean = 0;
	int itt = 0;
	int cnt = 1;

	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());
	rng.seed(RUN_SEED_SAA[siter]);
	int eval_num = 0;
	double UBsum = 0;
	double UB2sum = 0;
	std::vector<IloNum> Cx(problem->stage_sub_prob[0].range_raw.getSize(), -99900999);

	while ((stdes > out_sample_std && itt <out_sample_max) || (itt<out_sample_min))
	{
		eval_num++;

		sampl_ing->set_random_rhs_EVAL(problem->stage_sub_prob[0].rhs, problem->SPprobINFO.RVs, rng);
		problem->solver.set_rhs_var(problem->stage_sub_prob[0], x,
			problem->stage_sub_prob[0].rhs, Cx);

		if (!problem->solver.Solve_Prob(problem->stage_sub_prob[0], true))
		{
			std::cout << "ERROR: SAA: EVAL: Prob Not Solved!" << std::endl;
			std::cout << problem->stage_sub_prob[0].model << std::endl;
			std::cout << "xFirst size: " << x.size() << std::endl;
			std::cout << "rhs size: " << problem->stage_sub_prob[0].rhs.size() << std::endl;
			for (int i = 0; i < x.size(); i++) std::cout << "master x" << i << ": " << x[i].value << std::endl;
			for (int i = 0; i < x.size(); i++)
				std::cout << "x" << i << ": " << problem->master_prob.sol[i].value << std::endl;
			std::cout << "master status: " << problem->master_prob.cplex.getStatus() << std::endl;
			std::cout << "sub status: " << problem->stage_sub_prob[0].cplex.getStatus() << std::endl;
			for (int i = 0; i < problem->stage_sub_prob[0].rhs.size(); i++)
				std::cout << "rhs" << i << ": " << problem->stage_sub_prob[0].rhs[i].value << std::endl;
			system("pause");
		}


		UB = problem->stage_sub_prob[0].zstar + problem->master_prob.zstar_wout_surro;

		UBsum += UB;
		UB2sum += UB*UB;
		
		if (isnan(UB))
		{
			std::cout << eval_num << " " << problem->stage_sub_prob[0].zstar << " , " << problem->master_prob.zstar_wout_surro << std::endl;
			std::cout << problem->stage_sub_prob[0].model << std::endl;
			std::cout << "UB: " << UB << ", UBsum: " << UBsum << ", UB2sum: " << UB2sum << std::endl;
			system("pause");
		}
		if (itt > 10) {
			UBmean = UBsum / (double)eval_num;
			stdes = std::sqrt(((double)1 / eval_num)*(abs((UB2sum / (double)eval_num) - UBmean * UBmean)));
		}
		itt++;

		if (itt < 500) {
			if (itt > cnt) {
				printf("-");
				cnt += 100;
			}
		}
		else {
			if (itt < 10000) {
				if (itt > cnt) {
					printf("+");
					cnt += 500;
				}
			}
			else {
				if (itt > cnt) {
					printf("#");
					cnt += 1000;
				}
			}
		}

	}

	rep_UB.push_back(UBmean);
	return itt;
}

void SAA::write_legend(std::ofstream& file, std::string probname)
{
	file << "--------------------------------------------------------------------------------" << std::endl;

	file << "                          SAMPLE AVERAGE APPROXIMATION                          " << std::endl;

	file << "--------------------------------------------------------------------------------" << std::endl;
	file << "--------------------------------------------------------------------------------" << std::endl;
	file << "--------------------------------------------------------------------------------" << std::endl;
	if (Regular_LShaped == 1) {
		file << "                 Regularized LShaped ALGORITHM IN EACH REPLICATION               " << std::endl;
	}
	else if (Regular_LShaped == 0) {
		file << "                        LShaped ALGORITHM IN EACH REPLICATION                    " << std::endl;
	}
	file << "--------------------------------------------------------------------------------" << std::endl;
	if (Sampling_Method == 0) {
		file << "Sampling: Monte Carlo" << std::endl;
	}
	else if (Sampling_Method == 1) {
		file << "Sampling: Latin Hypercube Sampling" << std::endl;
	}
	file << "--------------------------------------------------------------------------------" << std::endl;
	if (Multi_LShaped == 0) {
		file << "Single Cut LShaped" << std::endl;
	}
	else if (Multi_LShaped == 1) {
		file << "Multi Cut LShaped" << std::endl;
	}
	file << "--------------------------------------------------------------------------------" << std::endl;
	file << "--------------------------------------------------------------------------------" << std::endl;
	file << "Number of Replications  = " << SAA_rep_num << std::endl;
	file << "Number of Initial Scenarios = " << N << std::endl;
	file << "Number of Scenarios of Evaluation = " << out_sample_max << std::endl;
	file << "Name of the test instance:    " << probname << std::endl;
	file << "--------------------------------------------------------------------------------" << std::endl;
	file << "--------------------------------------------------------------------------------" << std::endl << std::endl;;
}

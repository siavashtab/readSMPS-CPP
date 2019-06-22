
#include"ProbPrep.h"


#pragma region Initialization_funcs
ProbPrep::ProbPrep()
{
	printHeader();
	SPprobINFO  = new SPProb_INFO();
	newSPparam();
	solver      = new Solver_CPLEX();
	master_prob = new Prob();
	mean_prob   = new Prob();
	newProb(*mean_prob);
}

ProbPrep::~ProbPrep()
{
	solver->end_solver();
	delete solver;
	delete master_prob;
	delete mean_prob;
	delete stage_sub_prob;
	delete SPprobINFO;
}

void ProbPrep::newRV(RV_info &rv)
{
	rv.CDF = new std::map<double, double>;
	rv.ColName = new std::string;
	rv.id = new std::string;
	rv.prob = new vec_d;
	rv.RowName = new std::string;
	rv.val = new vec_d;
}

void ProbPrep::newSPparam()
{
	SPprobINFO->RVs = new std::vector<RV_info>;
	SPprobINFO->STOC_TYPE = new std::string;
	SPprobINFO->TIME_col_idx = new vec_i;
	SPprobINFO->TIME_info = new std::vector<std::vector<std::string>>;
	SPprobINFO->TIME_row_idx = new vec_i;
}

void ProbPrep::newProb(Prob &tmpProb)
{
	tmpProb.beta = new vec2_d;
	tmpProb.beta_sum = new vec_d;
	tmpProb.Cx = new vec_d;
	tmpProb.duals = new std::vector<SOL_str>;
	tmpProb.dual_R = new vec_d;
	tmpProb.dual_tot = new vec_d;
	tmpProb.strType = new std::string;
	//tmpProb.feas_cut_duals = new std::vector<SOL_str>;
	tmpProb.id = new std::string;
	//tmpProb.LShaped_feas = new IloRangeArray;
	//tmpProb.master_cuts = new IloRangeArray;
	tmpProb.obj_coef_raw = new std::vector<Coeff_Sparse>;
	//tmpProb.opt_cut_duals = new std::vector<SOL_str>;
	//tmpProb.other_cuts = new std::vector<IloRange>;
	tmpProb.prev_rng_coefs_raw = new std::vector<std::vector<Coeff_Sparse>>;
	tmpProb.prev_vars_raw = new std::vector<IloNumVarArray>;
	tmpProb.probType = new std::string;
	tmpProb.rho = new vec_d;
	tmpProb.rhs = new std::vector<SOL_str>;
	tmpProb.rng_coefs_raw = new std::vector<std::vector<Coeff_Sparse>>;
	//tmpProb.R_vars_raw = new IloNumVarArray;
	tmpProb.r_w = new vec_d;
	tmpProb.sol = new std::vector<SOL_str>;
	tmpProb.sol_R = new std::vector<SOL_str>;
	tmpProb.sol_surrogate = new std::vector<SOL_str>;
	tmpProb.strType = new std::string;
	tmpProb.surro_idx = new vec2_i;
	tmpProb.xhat = new std::vector<SOL_str>;
}

void ProbPrep::initialize(std::string filename, int scen_num)
{
	
	SPprobINFO->file_name = filename;
	SPprobINFO->initial_scen_num = &scen_num;
	
	if (Algorithm == 0)
	{
		SPprobINFO->algorithm = "LShaped";
	}
	else
	{
		SPprobINFO->algorithm = "PH";
	}
	SPprobINFO->dirr_model = ".\\model\\" + SPprobINFO->file_name + "\\";
	SPprobINFO->dirr_algo = ".\\Algorithm\\" + SPprobINFO->algorithm + "\\" + SPprobINFO->file_name + "\\";
	SPprobINFO->dirr_input = ".\\Input\\" + SPprobINFO->file_name + "\\";
	SPprobINFO->dirr_saa = ".\\SAA\\" + SPprobINFO->algorithm + "\\" + SPprobINFO->file_name + "\\";
	SPprobINFO->algo_type = ".txt";
	SPprobINFO->model_type = ".lp";
	SPprobINFO->output_algo = SPprobINFO->file_name + "_" + std::to_string(*SPprobINFO->initial_scen_num) + "_";
	SPprobINFO->output_model = "model_" + SPprobINFO->file_name + "_";
	SPprobINFO->output_saa = "SAA_" + SPprobINFO->file_name + "_";

	if (Regular_LShaped == 1) {
		SPprobINFO->output_algo += "reg_";
	}

	if (Multi_LShaped == 1) {
		SPprobINFO->output_algo += "mult";
	}

	if (Regular_LShaped == 1) {
		SPprobINFO->output_saa += "reg_";
	}
	if (Sampling_Method == 0) {
		SPprobINFO->output_saa += "mont_";
	}
	else if (Sampling_Method == 1) {
		SPprobINFO->output_saa += "lhs_";
	}
	if (Multi_LShaped == 1) {
		SPprobINFO->output_saa += "mult_";
	}
	
	Probs_from_SMPS();
	
	if (Algorithm == 0)
	{
		*master_prob = create_master_prob();
		
		for (int st = 1; st < *SPprobINFO->stage_num; st++)
		{
			stage_sub_prob->push_back(create_sub_probs(st));
			//cout << " RHS size" << stage_sub_prob[0].rhs.size();
			//cout << stage_sub_prob[0].model << endl;
		}
	}//if algo is LShaped based
	else if (Algorithm == 1)
	{
		std::cout << "This code is for Lshaped Method" << std::endl;
	}//if algo is ph based



}

void ProbPrep::printHeader()
{
	if (SMPS_Prob == 1)
	{
		std::cout << "-------------------------------------------------------------------------------" << std::endl;
		std::cout << "              SP problem is going to be read from SMPS files                   " << std::endl;
		std::cout << "-------------------------------------------------------------------------------" << std::endl << std::endl;
	} //if reading from SMPS files
	else
	{
		std::cout << "-------------------------------------------------------------------------------" << std::endl;
		std::cout << "                           SP problem is going to be created                   " << std::endl;
		std::cout << "-------------------------------------------------------------------------------" << std::endl << std::endl;
	} //If probs written manually
}
#pragma endregion class initialization functions 


#pragma region SMPS_reader_section
//reading the file from SMPS files
void ProbPrep::Probs_from_SMPS()
{
	read_COR();
	read_TIM();
	read_STOC();
}

void ProbPrep::read_COR()
{
	solver->open_solver();
	
	std::string CORfileName = SPprobINFO->dirr_input + SPprobINFO->file_name + ".cor";
	std::cout << "looking for: " << CORfileName << std::endl;
	mean_prob->name = CORfileName.c_str();
	*mean_prob->strType = "mean";
	solver->open_prob(*mean_prob);
	solver->ImportModel(*mean_prob);
	
	if (false) printf(" ** flag ** \n");
	if (show_meanprob == 1) std::cout << mean_prob->model << std::endl;
	
	solver->get_Linear_obj_coeffs(*mean_prob);
	solver->get_Linear_rng_coeffs(*mean_prob);
	
	SPprobINFO->num_var = mean_prob->num_var;
	SPprobINFO->num_rngs = mean_prob->num_rng;

	mean_prob->opt = false;
	mean_prob->feas = false;
	
	mean_prob->isReg = false;
	solver->Create_Prob(*mean_prob);
	solver->Solve_Prob(*mean_prob, false);
	
	
	mean_prob->model_name = SPprobINFO->dirr_model + SPprobINFO->output_model +
		"mean" + SPprobINFO->model_type;
	mean_prob->name = mean_prob->model_name.c_str();

	if (print_lp_models == 1) solver->Print_Prob(*mean_prob);

	std::cout << "Objective Function of the Mean Problem:  " << mean_prob->zstar << std::endl;
}

void ProbPrep::read_TIM()
{
	std::string TIMfileName = SPprobINFO->dirr_input + SPprobINFO->file_name + ".tim";
	std::ifstream TIMfileStream;
	TIMfileStream.open(TIMfileName.c_str());

	if (TIMfileStream.fail())
	{
		std::cout << "                 ERROR: Could Not Open TIM file" << std::endl;
		exit(1);
	}

	std::string line;
	int stage = 0;
	while (getline(TIMfileStream, line)) {
		if (line[0] == '*' || line[0] == '\r' || (line[0] >= '0' && line[0] <= 'Z'))
			continue; 

		std::stringstream lineStream(line);
		std::string col_name, row_name;
		lineStream >> col_name >> row_name;		

		SPprobINFO->TIME_info->push_back(std::vector<std::string>(2));

		SPprobINFO->TIME_info->at(stage)[0] = col_name;
		SPprobINFO->TIME_info->at(stage)[1] = row_name;

		stage++;							  
	}

	SPprobINFO->stage_num = &stage;

	TIMfileStream.close();
}

//read STOC file
void ProbPrep::read_STOC()
{
	std::string stoFileName = SPprobINFO->dirr_input + SPprobINFO->file_name + ".sto";
	std::ifstream stoFileStream;
	stoFileStream.open(stoFileName.c_str());
	if (stoFileStream.fail()) {
		std::cout << "                 ERROR: Could Not Open STOC file" << stoFileName << std::endl;
		exit(1);
	}

	std::string col1, col2, col3, col4;
	int flag = 0;

	while (!safeGetline(stoFileStream, col1).eof())
	{

		std::stringstream lineStream(col1);
		lineStream >> col1 >> col2 >> col3 >> col4;

		if (col1.compare("INDEP") == 0) // getting new scenarios
		{
			std::cout << "                 READING INDEP PROP DISTRIBUTION" << std::endl;
			break;
		}

	}

	if (flag == 0)
	{
		stoFileStream.close();
		read_STOC_INDEP();
	}//the distribution is indep discrete (use read_STOC_INDEP())


}

//read STOC file
void ProbPrep::read_STOC_INDEP()
{
	std::string stoFileName = SPprobINFO->dirr_input + SPprobINFO->file_name + ".sto";
	std::ifstream stoFileStream;
	stoFileStream.open(stoFileName.c_str());

	std::string col1, col2, col3, col4;
	int flag = 0;
	RV_info rv_empty;
	newRV(rv_empty);
	double sum = 0;
	
	while (!safeGetline(stoFileStream, col1).eof())
	{

		std::stringstream lineStream(col1);
		lineStream >> col1 >> col2 >> col3 >> col4;
		
		if (col1.compare("INDEP") == 0) // getting new scenarios
		{
			*SPprobINFO->STOC_TYPE = "INDEPENDENT DISCRETE";
		}
		else if (col3.compare("") != 0)
		{
			if (col2.compare("DISCRETE") == 0)
			{
				flag++;
			}
			else if (flag == 1)
			{
				*rv_empty.ColName = col1;
				*rv_empty.RowName = col2;
				rv_empty.val->push_back(atof(col3.c_str()));
				rv_empty.prob->push_back(atof(col4.c_str()));
				sum += atof(col4.c_str());
				rv_empty.CDF->insert(std::make_pair(atof(col3.c_str()), sum));
				flag = INFINITY;
			}
			else if (col1.compare("*") == 0 || col1.compare("ENDATA") == 0)
			{
				SPprobINFO->rv_num++;
				SPprobINFO->RVs->push_back(rv_empty);
				rv_empty.val->clear();
				rv_empty.prob->clear();
				rv_empty.CDF->clear();
				sum = 0;
			}
			else
			{
				*rv_empty.ColName = col1;
				*rv_empty.RowName = col2;
				rv_empty.val->push_back(atof(col3.c_str()));
				rv_empty.prob->push_back(atof(col4.c_str()));
				sum += atof(col4.c_str());
				rv_empty.CDF->insert(std::make_pair(atof(col3.c_str()), sum));
			}
		}


		if (col1.compare("ENDATA") != 0) continue;

		for (int rv = 0; rv < SPprobINFO->RVs->size(); rv++)
		{
			if (*SPprobINFO->RVs->at(rv).ColName == "RHS")
			{
				SPprobINFO->RVs->at(rv).rv_col_num = 0;
				for (int rg = 0; rg < solver->getSize(*mean_prob->range_raw); rg++)
				{
					if ((std::string)solver->getName(mean_prob->range_raw->operator[](rg)) == *SPprobINFO->RVs->at(rv).RowName)
					{
						SPprobINFO->RVs->at(rv).rv_row_num = rg;
						SPprobINFO->RVs->at(rv).name = solver->getName(mean_prob->range_raw->operator[](rg));
					}//if range name is equal
				}//for ranges
			}//rv in RHS
			else
			{
				printf("              ERROR: RV is not in RHS    \n");
				std::cout << *SPprobINFO->RVs->at(rv).ColName << std::endl;
			}//rv is not in RHS
		}//for all rvs

	}

}

#pragma endregion reading and storing the info form SMPS files



#pragma region Master_Problem_Creation
Prob ProbPrep::create_master_prob()
{
	Prob prob;
	newProb(prob);

	solver->open_prob(prob);
	
	add_surrogate_master_vars(prob, *solver);  //Add surrogate LShaped vars for defining lower bounding cuts
	add_master_vars(prob, *solver);            //separate master variables from mean prob and add to master_prob
	
	*prob.strType = "master";
	
	for (int v1 = 0; v1 < mean_prob->obj_coef_raw->size(); v1++)
	{
		for (int v2 = 0; v2 < solver->getSize(*prob.vars_raw); v2++)
		{
			if (mean_prob->obj_coef_raw->at(v1).col_name == solver->getName(prob.vars_raw->operator[](v2)))
			{
				Coeff_Sparse empt;
				empt = mean_prob->obj_coef_raw->at(v1);
				empt.col = v2;
				prob.obj_coef_raw->push_back(empt);
				break;
			}
		}
	}
	
	if (Regular_LShaped == 1)
	{
		prob.lambda = 0.0;
		prob.sigma = LShaped_reg_sigma;
	}
	else
	{
		prob.lambda = 0.0;
		prob.sigma = 0.0;
		prob.isReg = false;
	}

	x_reg.resize(prob.vars_raw->getSize());
	for (int i = 0; i < prob.vars_raw->getSize(); i++) x_reg[i].value = mean_prob->sol->at(i).value;

	add_master_rngs(prob, *solver);
	add_master_obj(prob, *solver);

	for (int i = SPprobINFO->TIME_row_idx->at(0); i < SPprobINFO->TIME_row_idx->at(1); i++)
	{
		prob.rng_coefs_raw->push_back(mean_prob->rng_coefs_raw->at(i));
	}
	
	prob.model_name = SPprobINFO->dirr_model + SPprobINFO->output_model + 
		                                    "master" + SPprobINFO->model_type;
	prob.name = prob.model_name.c_str();

	prob.env->setDeleter(IloLinearDeleterMode);

	solver->Create_Prob(prob);

	if (print_lp_models == 1) print_master_prob(prob);
	int s1 = solver->getSize(*prob.range_raw);
	int s2 = solver->getSize(*prob.vars_raw);
	prob.num_rng = &s1;
	prob.num_var = &s2;
	*prob.probType = solver->getProbtype(prob);
	prob.random = false;
	if (Regular_LShaped == 1) prob.isReg = true;

	solver->Solve_Prob(prob, true);
	std::cout << "Empty Master Objective:  " << prob.zstar << std::endl;
	std::cout << "Master Type:  " << solver->getProbtype(prob) << std::endl;


	//Master problem does not have any cuts at this point
	prob.opt = false;
	prob.feas = false;
	prob.update_xstar = false;

	std::cout << " ----- SUCCESS: master is completed ----- " << std::endl;
	std::cout << " ---------------------------------------- " << std::endl;

	return prob;
}

void ProbPrep::add_surrogate_master_vars(Prob& prob, Solver_CPLEX& solver)
{
	prob.has_surrogate = true;
	
	if (prob.surro_vars_raw->getSize() > 0)  prob.surro_vars_raw->clear();
	prob.surro_vars_raw = new IloNumVarArray(*prob.env);

	//first add surrogate variables
	for (int s = 0; s < *SPprobINFO->initial_scen_num; s++)
	{
		char varname[20];
		sprintf(varname, "\eta_%d", (int)(s+1));
		solver.add_to_array(*prob.surro_vars_raw, 
			              solver.Create_Var_explicit(-IloInfinity, IloInfinity, 1, varname));
	}

	
	std::cout << " ----- Surrogate vars are added ----- " << std::endl;
}

void ProbPrep::add_master_vars(Prob& prob, Solver_CPLEX& solver)
{
	int stage_count = 0;
	for (int v = 0; v < mean_prob->vars_raw->getSize(); v++)
	{
		if (SPprobINFO->TIME_info->at(0)[0] == (std::string)solver.getName(mean_prob->vars_raw->operator[](v)))
		{
			SPprobINFO->TIME_col_idx->push_back(v);
			solver.add_to_array(*prob.vars_raw, mean_prob->vars_raw->operator[](v));
			stage_count++;
		}
		else if (SPprobINFO->TIME_info->at(stage_count)[0] == (std::string)solver.getName(mean_prob->vars_raw->operator[](v)))
		{
			SPprobINFO->TIME_col_idx->push_back(v);
			if (SPprobINFO->TIME_info->size() == stage_count + 1)
			{
				break;
			}
			else
			{
				stage_count++;
			}		
		}
		else
		{
			solver.add_to_array(*prob.vars_raw, mean_prob->vars_raw->operator[](v));
		}
	}
	
	std::cout << " ----- master vars are added ----- " << std::endl;
}

void ProbPrep::add_master_rngs(Prob& prob, Solver_CPLEX& solver)
{
	int stage_count = 0;
	if (SPprobINFO->TIME_info->at(0)[1] == (std::string)solver.getName(*mean_prob->obj_raw))
	{
		stage_count++;
		SPprobINFO->TIME_row_idx->push_back(0);
	}

	for (int r = 0; r < solver.getSize(*mean_prob->range_raw); r++)
	{
		if (SPprobINFO->TIME_info->at(0)[1] == (std::string)solver.getName(mean_prob->range_raw->operator[](r)))
		{
			SPprobINFO->TIME_row_idx->push_back(r);
			solver.add_to_array(*prob.range_raw, mean_prob->range_raw->operator[](r));
			stage_count++;
		}
		else if(SPprobINFO->TIME_info->at(stage_count)[1] == (std::string)solver.getName(mean_prob->range_raw->operator[](r)))
		{
			SPprobINFO->TIME_row_idx->push_back(r);
			if (SPprobINFO->TIME_info->size() == stage_count + 1)
			{
				break;
			}
			else
			{
				stage_count++;
			}
		}
		else
		{
			solver.add_to_array(*prob.range_raw, mean_prob->range_raw->operator[](r));
		}
	}
	SPprobINFO->TIME_row_idx->push_back(solver.getSize(*mean_prob->range_raw));
	//cout << " ----- master rngs are added ----- " << endl;
}

void ProbPrep::add_master_obj(Prob& prob, Solver_CPLEX& solver)
{
	IloExpr empty_expr(*prob.env);

	empty_expr = prob.surro_vars_raw->operator[](0);

	if (Multi_LShaped == 1)
	{
		for (int s = 1; s < prob.surro_vars_raw->getSize(); s++)
		{
			empty_expr += prob.surro_vars_raw->operator[](s);

		}
	}

	double minVal = 1/(double)*SPprobINFO->initial_scen_num;

	if (prob.obj_coef_raw->size() > 0)
	{
		for (int v = 0; v < prob.obj_coef_raw->size(); v++)
		{
			empty_expr += prob.obj_coef_raw->at(v).val * prob.vars_raw->operator[](prob.obj_coef_raw->at(v).col);
			if (prob.obj_coef_raw->at(v).val < minVal) minVal = prob.obj_coef_raw->at(v).val;
		}
	}
	
	//prob.model.add(empty_expr >= mean_prob.zstar);
	if (Regular_LShaped == 1)
	{
		if (minVal == 1 / (double)*SPprobINFO->initial_scen_num)
		{
			prob.lambda = abs(minVal) / 50;
			prob.sigma  = abs(minVal) / 100;
		}
		else
		{
			prob.lambda = abs(minVal) / 2;
			prob.sigma = abs(minVal) / 4;
		}


		empty_expr += solver.set_QP_obj(prob, x_reg);
	}
	

	*prob.obj_raw = solver.Create_OBJ_explicit(mean_prob->obj_raw->getSense(), empty_expr, mean_prob->obj_raw->getName());
	empty_expr.end();
	//cout << " ----- master obj is added ----- " << endl;
}

void ProbPrep::print_master_prob(Prob& prob)
{
	solver->Print_Prob(prob);
	//cout << " ----- master is printed ----- " << endl;
}
#pragma endregion creating the master problem


#pragma region stage_based_subprob_generation
Prob ProbPrep::create_sub_probs(int stage)
{
	Prob prob;
	newProb(prob);

	solver->open_prob(prob);
	*prob.strType = "sub";
	
	int end;
	if (SPprobINFO->TIME_col_idx->size() > stage + 1)
	{
		end = SPprobINFO->TIME_col_idx->at(stage + 1);
	}
	else
	{
		end = solver->getSize(*mean_prob->vars_raw);
	}

	add_sub_vars(prob, *solver, SPprobINFO->TIME_col_idx->at(stage), end);
	std::cout << " ----- sub_"<< stage <<" vars are added ----- " << std::endl;
	add_sub_obj(prob, *solver);
	std::cout << " ----- sub_" << stage << " obj is added ----- " << std::endl;
	add_sub_rngs(prob, *solver);
	std::cout << " ----- sub_" << stage << " ranges are added ----- " << std::endl;

	prob.model_name = SPprobINFO->dirr_model + SPprobINFO->output_model +
		"sub" + std::to_string(stage) + SPprobINFO->model_type;

	solver->Create_Prob(prob);
	prob.isReg = false;

	if (print_lp_models == 1) print_sub_prob(prob, stage);
	//cout << " ----- sub_" << stage << " is printed ----- " << endl;
	int n1 = solver->getSize(*prob.range_raw);
	int n2 = solver->getSize(*prob.vars_raw);
	prob.num_rng = &n1;
	prob.num_var = &n2;
	*prob.probType = solver->getProbtype(prob);
	prob.random = true;

	fill_out_rhs(prob, *solver);

	prob.opt = false;
	prob.feas = false;

	std::cout << " ----- SUCCESS: sub_" << stage <<" is completed ----- " << std::endl;
	std::cout << " ---------------------------------------- " << std::endl;

	return prob;
}

void ProbPrep::add_sub_vars(Prob& prob, Solver_CPLEX& solver, int start, int end)
{
	for (int i = start; i < end; i++)
	{
		solver.add_to_array(*prob.vars_raw, mean_prob->vars_raw->operator[](i));
	}
}

void ProbPrep::add_sub_rngs(Prob& prob, Solver_CPLEX& solver)
{

	solver.decompose_range(*mean_prob->range_raw, prob, *prob.vars_raw, 
		                   SPprobINFO->TIME_row_idx->at(1), SPprobINFO->TIME_row_idx->at(2));
}

void ProbPrep::add_sub_obj(Prob& prob, Solver_CPLEX& solver)
{
	
	for (int v = 0; v < mean_prob->obj_coef_raw->size(); v++)
	{
		for (int v2 = 0; v2 < solver.getSize(*prob.vars_raw); v2++)
		{
			if (mean_prob->obj_coef_raw->at(v).col_name == solver.getName(prob.vars_raw->operator[](v2)))
			{
				Coeff_Sparse empt;
				empt = mean_prob->obj_coef_raw->at(v);
				empt.col = v2;
				prob.obj_coef_raw->push_back(empt);
				break;
			}
		}
	}
	
	IloExpr empt_expr(*prob.env);

	if (prob.obj_coef_raw->size() > 0)
	{
		for (int i = 0; i < prob.obj_coef_raw->size(); i++)
		{
			empt_expr += prob.obj_coef_raw->at(i).val * prob.vars_raw->operator[](prob.obj_coef_raw->at(i).col);
		}
	}

	*prob.obj_raw = solver.Create_OBJ_explicit(mean_prob->obj_raw->getSense(), empt_expr, mean_prob->obj_raw->getName());

}

void ProbPrep::print_sub_prob(Prob& prob, int stage)
{
	solver->Print_Prob(prob);
}

void ProbPrep::fill_out_rhs(Prob& prob, Solver_CPLEX& solver)
{
	for (int r = 0; r < *prob.num_rng; r++)
	{
		SOL_str empt;
		empt.row = solver.getName(prob.range_raw->operator[](r));
		empt.row_num = r;
		empt.value = solver.get_rhs(prob.range_raw->operator[](r));
		empt.isRandom = false;
		empt.randomVar_indx = -1;
		prob.rhs->push_back(empt);
	}

	for (int rv = 0; rv < SPprobINFO->RVs->size(); rv++)
	{
		for (int r = 0; r < *prob.num_rng; r++)
		{
			if (*SPprobINFO->RVs->at(rv).RowName == (std::string)prob.rhs->at(r).row)
			{
				prob.rhs->at(r).isRandom = true;
				prob.rhs->at(r).randomVar_indx = rv;
				break;
			}
		}
	}

	if (sub_rhs_print == 1)
	{
		for (int r = 0; r < *prob.num_rng; r++)
		{
			std::cout << prob.rhs->at(r).row << " ";
			std::cout << prob.rhs->at(r).row_num << " ";
			std::cout << prob.rhs->at(r).isRandom << " ";
			std::cout << prob.rhs->at(r).randomVar_indx << " ";
			std::cout << prob.rhs->at(r).value << std::endl;
		}
	}
}
#pragma endregion creating stage based decomposed subproblems


#pragma region PH_subprob_creation

#pragma endregion creating subproblems needed for PH algorithm
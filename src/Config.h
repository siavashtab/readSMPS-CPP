#ifndef Config_H
#define Config_H

//RUNNING PARAMS
#define   EPGAP         0.0001
#define   TimeLim        100
#define   Parallel      0                       // (1) for parallel computing
#define   EPs           1e-10
#define   num_stages    2                       // number of stages (2 for two stage SP)
#define   OPt_GAP       0.01                    // optimality gap
#define   numDigits     2                       // number of decimals



//SAA Infos
#define   Experiments     1                      //doing experiments on number of scenarios?
#define   Expr_num        5                      //number of experiments
#define   SAA_rep_num     30                     //SAA number of replication
#define   initial_scen    10                    //initial number of scenarios
#define   out_sample_std    0.002                //out sample error
#define   out_sample_max    20000                //out sample error
#define   out_sample_min    1000                 //out sample error


//Printing options
#define   sub_rhs_print          1      //for printing the information of the subproblem rhs
#define   print_lp_models        1      //1(print) - 0(not print)
#define   show_meanprob          0      //display the mean problem
#define   show_model_scen        0      //display the model of each scenarios
#define   show_model_scen_a      0      //display the model of each scenarios(Adaptive)
#define   show_master_iter       0      //display master model after each iteration
#define   print_Algo_results     1      //pront the algorthm results
#define   print_Master_At_end     0     //display the master problem after LShaped with all the cuts
#define   display_LShaped_LB_UB       0     
#define   display_estimator_iter0     0
#define   display_Master_surrogates	  0



//Sampling Configs
#define   Sampling_Method                1           //monte carlo(0) - lhs(1) - ACS(2)
#define   Display_Samples                0           //display the taken samples

//File Reading Configs
#define   SMPS_Prob          1      //smps format prob(1) or a written prob(0)

//Algorithm Configs
#define   Algorithm              0      //LShaped(0) - PH(1)
#define   Regular_LShaped        1      //Regularized LShaped(1) - Linear obj master LShaped (1)
#define   Multi_LShaped          1      //Single(0) - Multi(1)
#define   Remove_Old_Cuts        0      //Remove old cuts
#define   LShaped_maxiter        500    
#define   LShaped_miniter        5 
#define   LShaped_gap            0.0001
#define   LShaped_reg_sigma      0.001

//ACS Configs
#define   ACS_cov           0.00001  //ACS neighborhood coverage in 0-1 interval
#define   min_n_c           20       //minimum # of samples within cluster
#define   max_n_c           1000     //maximum # of samples within cluster pgp2(50) storm(50) ssn(50) retail(100)
#define   neighb            1000     //fraction of answers in the neighborhood strom(0.001) pgp2(1) ssn(0.01) retail(10)
#define   min_thre          1        //minimum # of samples within cluster
#define   disp_point_estim           1        //display the point estimate
#define   disp_point_estim_iter      0        //display the point estimate
#define   print_point_Estim_samples  0  
#define   print_SSW_info             0
#define   print_SSW_info_seq         0
#define   scen_or_clust              0        //scenario based(0) ACS or cluster based(1)
#define   scen_type                  1        //monte(0) lhs(1)


#endif;

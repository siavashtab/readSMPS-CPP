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



//Sampling Configs
#define   Sampling_Method                1           //monte carlo(0) - lhs(1) - ACS(2)
#define   Display_Samples                0           //display the taken samples

//File Reading Configs
#define   SMPS_Prob          1      //smps format prob(1) or a written prob(0)

//Algorithm Configs
#define   Algorithm              0      //LShaped(0) - PH(1)
#define   Regular_LShaped        1      //Regularized LShaped(1) - Linear obj master LShaped (1)
#define   Multi_LShaped          1      //Single(0) - Multi(1)



#endif;

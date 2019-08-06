# Stochastic Control & Verification via Barrier Functions
This repository contains the code for the case studies in the paper "Verification  and  Control  for  Finite-Time  Safety  ofStochastic  Systems  via  Barrier  Functions" (CCTA 2019) and "A Barrier Function Approach to Finite-Time StochasticSystem Verification and Control" _(in submission)_

Author: Cesar Santoyo <br />
E-mail: csantoyo@gatech.edu

If you have any questions, please e-mail the author at the e-mail above.

## Required Software: ##
* MATLAB
* SOSTOOLS (https://www.cds.caltech.edu/sostools/)
* SDPT3 (http://www.math.nus.edu.sg/~mattohkc/sdpt3.html)

Note: These case studes were constructed using MATLAB 2018. Additionally, the symbolic toolbox is required. It is possible to execute the code without the symbolic toolbox (see SOSTOOLs manual for details).

## Case Study 1: ##
cs1_main.m:  <br />
* Running this file will run the algorith which was used to produce the results of this paper. The required toolboxes are mentioned above. You may run the individual dependencies separately to get a more careful look at the corresponding results. This case study is for the 1-D stochastic dynamics. 

## Case Study 2: ##
cs2_main.m:  <br />
* Running this file will run the algorithm which was used to produce the results of this paper. The required toolboxes are mentioned above. You may run the individual dependencies separately to get a more careful look at the corresponding results. These results are for the 2-D stochastic dynamics. 

## Case Study 3: ##
cse_main.m:  <br />
* Running this file will run the algorithm which was used to produce the results of the discrete-time example in the journal paper _(in submission)_. The required toolboxes are mentioned above. You may run the individual dependencies separately to get a more careful look at the corresponding results. These results are for the 2 state population growth model. 

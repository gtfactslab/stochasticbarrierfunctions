%           Main file for generating plots for the third case study in the 
%           Automatica Paper (in submission). 
%               
% Paper:    A BarrierFunction Approach to Finite-Time Stochastic System Verification and Control
%
% Author:   Cesar G. Santoyo
% Date:     August 6th, 2019
% 
% If you have any questions, e-mail the author at csantoyo@gatech.edu
%
% Description:  This script runs the results demonstrated in the Automatica
%               paper for Caset Study 3.


%% Run verification algorithm (Algo 1 from paper)
run('cs3algo1_TAC.m');


%% Run control synthesis algorithm
% This portion is for control synthesis using a linear barrier function in
% the positive domain
%
% For this uncomment lines: 27, 28, 31, 32, 37, 41, 43
% For this comment lines: 36, 40, 42
%
run('cs3algo3_TAC');
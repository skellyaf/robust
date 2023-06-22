close all; clear; clc;

addpath(genpath('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm'));
addpath('C:\Users\skell\OneDrive - USU\Documents\code_repos\robust\Analysis\internal_utilities');
cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221204_0956.31_leoCrit_to_leoSunSync_robustMissionSpecific_arrivalJ')
% cd('C:\Users\skell\OneDrive - USU\Documents\code_repos\multiseg_opt_manual_tcm\sims\20221207_1450.18_leo2meo_robust_fixed_tf')

load('workspace.mat')
close all;

simparams.dynSys = '2bp';

plotIterationHistory(history.x, simparams);


saveallfigs('./',0);
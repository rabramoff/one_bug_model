close all;
clear all;
clc;

%remark:
%in all those experiments, Ms is set to 5 gC equivalent. 
[status,results]=system('pwd');
sstrs=strsplit(results,'/one_bug_model');
addpath(genpath([sstrs{1},'/one_bug_model/']));
Is=0.00015*24*365; %SOC input gC/year
Ic=0.00001*24*365; %DOC input gC/year
%one_box_deb_transtT_Input_driver(Is,Ic,Is,0);

one_box_deb_transtT_Input_driver(Is,Ic,15*Is,15*Ic);

one_box_deb_transtT_Input_driver(Is,Ic,20*Is,20*Ic);
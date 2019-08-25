clear all
clc
close all
global Ph FT Dp v m_pl P_pl Cstr C MP cx Shape Mat DesignParameters DsPm MPr ...
    Motor_forest DOE_forest Mdl_pr_forest DOE_pr_forest prop_pr_forest...
    prop_pr_DOE   err_flag output t count INPUT RUN results FUNoutput input

Missions=["chain1","chain2","chain3","RussianDoll3","Relay_Mother"];

count=1;
RUN=1;


for i=1:length(Missions)
    [~]=level1(Missions(i));
end

results;  
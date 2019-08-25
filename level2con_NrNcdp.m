%%% This function is a wrap of the LaWorkout function, it outputs the
%%% constrainsts.

function [cons,ceq]=level2con_NrNcdp(input)
 global INPUT
% % define the design parameters to be used in sub-level functions
%                     DsPm(1)=input(1)*2;       %N_r(i);
%                     DsPm(2)=input(2);       %N_c(j);
%                     DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
%                     DsPm(4)=input(4)*0.0254;       %d_pm;
%                     DsPm(5)=input(3);       %Mat_body(l);
%                     DsPm(6)=input(1)*2;       % N_arm
%                     DsPm(7)=input(6);       % TrW
%                     DL=input(5);
input1=[INPUT,input];
FUNoutput=level2_opt(input1);
err=FUNoutput(19);

cons=err-0.05;      % if the error can't converge to 0, at least,
% the error should not be larger than 2% of the initial GTOM guess. 

ceq=[];

end
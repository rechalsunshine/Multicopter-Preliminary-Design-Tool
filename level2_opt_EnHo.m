%%% This function is a wrap of the level2_opt function, it outputs the MOTM

function EnHo=level2_opt_EnHo(input)
global  t result
% % define the design parameters to be used in sub-level functions
%                     DsPm(1)=input(1)*2;       %N_r(i);
%                     DsPm(2)=input(2);       %N_c(j);
%                     DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
%                     DsPm(4)=input(4)*0.0254;       %d_pm;
%                     DsPm(5)=input(3);       %Mat_body(l);
%                     DsPm(6)=input(1)*2;       % N_arm
%                     DsPm(7)=input(6);       % TrW
%                     DL=input(5);

FUNoutput=level2_opt(input);
if abs(FUNoutput(19))<0.02
%     MTOM=output(18);
% price=FUNoutput(20);
% Energy=output(21);
EnHo=FUNoutput(23);
% EnCru=output(24);
% Vfw=output(25);
% Vas=output(26);
result(t,:)=FUNoutput;
t=t+1;
fprintf('current hovering endurance: %f min, error: %.1f%% \n', EnHo, FUNoutput(19)*100);
else
%         MTOM=output(18)*1000;
% price=FUNoutput(20)*1000;
% Energy=output(21)*1000;
EnHo=FUNoutput(23)*0.01;
% EnCru=output(24)*1000;
% Vfw=output(25)*1000;
% Vas=output(26)*1000;
end

EnHo=100-EnHo;
end
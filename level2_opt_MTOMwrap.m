%%% This function is a wrap of the level2_opt function, it outputs the MOTM

function MTOM=level2_opt_MTOMwrap(input)
global  DsPm FUNoutput result t
% define the design parameters to be used in sub-level functions
                    DsPm(1)=round(input(1)/2)*2;       %N_r(i);
                    DsPm(2)=input(2);       %N_c(j);
                    DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
                    DsPm(4)=input(4)*0.0254;       %d_pm;
                    DsPm(5)=input(3);       %Mat_body(l);
                    DsPm(6)=input(1)*2;       % N_arm
                    DL=input(5);
                    err=level2_NrNcdp(DL);
results=FUNoutput;
if err<0.05
    MTOM=results(18);
% price=output(20);
% Energy=output(21);
% EnHo=output(23);
% EnCru=output(24);
% Vfw=output(25);
% Vas=output(26);
result(t,:)=results;
t=t+1;
err_100=results(19)*100;
fprintf('current MTOM: %f g, error: %.1f%% \n', MTOM, err_100);
else
        MTOM=FUNoutput(18)*1000;
        fprintf('invalid results...\n');
% price=output(20)*1000;
% Energy=output(21)*1000;
% EnHo=output(23)*1000;
% EnCru=output(24)*1000;
% Vfw=output(25)*1000;
% Vas=output(26)*1000;
end
end
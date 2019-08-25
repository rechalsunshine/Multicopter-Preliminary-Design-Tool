%%% This function is a wrap of the level2_opt function, it outputs the MOTM

function MTOM=level2_opt_MTOM(input)
global  result t DsPm Cstr FUNoutput
% define the design parameters to be used in sub-level functions
                    DsPm(1)=round(input(1)/2)*2;       %N_r(i);
                    DsPm(2)=input(2);       %N_c(j);
                    DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
                    DsPm(4)=input(4)*0.0254;       %d_pm;
                    DsPm(5)=input(3);       %Mat_body(l);
                    DsPm(6)=DsPm(1);       % N_arm

                    
                    fun_dl=@level2_NrNcdp;    % set the optimisation function
% Disk Loading is the only variable to be changed
ub_dl=Cstr(17);    % set upper boundary
lb_dl=Cstr(16);    % set lower boundary
DL0=lb_dl+0.001;  % set the starting point
options=optimset('Display','iter','TolFun',1e-6,'TolX',1e-5);
[DL_opt,Min_err_per,existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
err=level2_NrNcdp(DL_opt);
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
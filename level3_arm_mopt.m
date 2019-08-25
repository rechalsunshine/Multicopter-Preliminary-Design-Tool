function [m_arm,pr_arm, opt_geo]=level3_arm_mopt(L_a)

global  Shape Mat input MPr

%%%% the strategy here is a conbined optimisation that uses GA to locate the local valley, and
% then use fmincon to go down the gradient. By doing so, it can find the accurate minimum and
% saves a lot of time (from 20 s to 0.3 s for a single run).

%%% constraints for ga optimisation

%%% constraints for fmincon optimisation
% Linear constraints, thickness does not exceed 1/2 of the overall size
% therefore, 2*property(4)-property(2)<=0, 2*property(4)-property(3)<=0


%%%%%  optimisation process %%%%
%1. GA for rough optimisation
optionsga=optimoptions('ga','PopulationSize',100,'MaxGenerations',300,...
    'MaxStallGenerations', 20,'FunctionTolerance',1e-5,'Display','off');
intcon=[1,6];                       % Integer variables
nonlconga=@level3_armcons;
lbga=[min(Shape),0.002,0.002,0.0008,L_a,min(Mat)];   % lower boundary of variables
ubga=[max(Shape),L_a/3,L_a/3,L_a/12,L_a,max(Mat)];

[geo_optga(1:6),ma_optga,exitflag] =...
    ga(@level3_armcalc_mid,6,[],[],[],[],lbga,ubga,nonlconga,intcon,optionsga);
if exitflag~=1
    optionsga=optimoptions('ga','PopulationSize',100,'MaxGenerations',500,...
    'MaxStallGenerations', 100,'FunctionTolerance',1e-5,'Display','off');
[geo_optga(1:6),ma_optga] =...
    ga(@level3_armcalc_mid,6,[],[],[],[],lbga,ubga,nonlconga,intcon,optionsga);

end
%2. fmincon for accurate optimisation
A=[-1,0,2,0; 0,-1,2,0];
b=[0; 0];
lb=[0.002,0.002,0.0008,L_a];                   % lower boundary of variables
ub=[L_a/5,L_a/5,L_a/12,L_a];        % ??? upper boundary of variables, relate the size to La_min
%x0=ub;
nonlcon=@level3_armcons_fmin;                   % nonlinear constraints
optionslocal=optimoptions('fmincon','Algorithm','sqp-legacy'  ,'CheckGradients', true,'OptimalityTolerance',...
    1e-8,'StepTolerance', 1e-8,'ConstraintTolerance',1e-7,'MaxFunctionEvaluations',...
    3000, 'MaxIterations',1000,'Display','off');      %'Diagnostics', 'on',,'Display','iter-detailed'
x0=geo_optga(2:5);    % use the result as starting point for fmincon
input(1)=geo_optga(1);
input(6)=geo_optga(6);
[geo_opt,ma_opt] = fmincon(@level3_armcalc_fmin,x0,A,b,[],[],lb,ub,nonlcon,optionslocal);
opt_geo=[input(1),geo_opt,input(6)];
m_arm=ma_opt;       % output mass
unit_price=MPr(opt_geo(1),opt_geo(6));
pr_arm=m_arm*unit_price;    % output price


%opt_geo(i*s,:)=[input(1),geo_opt(i*s,:),input(6)];
%t=11;



%M_arm=min(ma_opt);

% Arm=LaWorkoutMax(opt_geo,load);
% if Arm.frequency_error_x<1e-7 && Arm.frequency_error_y<1e-7 &&...
%         Arm.deflection_error<1.2e-7 && Arm.shear_stress<1e-7 &&...
%         Arm.normal_compression_stress<1e-7 && Arm.normal_tension_stress<1e-7
%     
%     disp("good optimisation");
% else
%     disp("!!!invalid result!!!");
%     %end
%     %runtime
% end
% 
% diag=2*La;          % diagonal distance equals twice of the minimum arm length, unit: m
% d_cen=diag*0.3172;           % central frame size is calculated as a certain ratio of the diagonal distance, unit: m
% M_cen=1.744*d_cen-0.3766/1000;     % central frame mass is correlated about the central frame size. unit: kg
% 
% M_lan=0.0293*MTOM0-1.82/100000;           %??? temperarily calculate the landing gear mass from the MTOM
% Me=M_arm*Narm+M_cen+M_lan;
% 
% 
% output(11+s)=d_cen*1000; % feature size of the central frame, unit; mm
% output(12+s)=M_cen*1000; % mass of central part, unit: g
% output(19+s)=M_arm*1000*Narm; % mass of one arm, unit: g
% output(22+s)=opt_geo(2)*1000; % arm cross section height, unit; mm
% output(23+s)=opt_geo(3)*1000; % arm cross section width, unit; mm
% output(24+s)=opt_geo(4)*1000; % arm cross section thickness, unit; mm
% output(26+s)=M_lan*1000; % mass of landing gears, unit: g

end

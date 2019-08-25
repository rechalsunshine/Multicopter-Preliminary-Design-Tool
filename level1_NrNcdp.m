clear all
clc
close all
global Ph FT Dp v m_pl P_pl Cstr C MP cx Shape Mat DesignParameters DsPm MPr ...
    Motor_forest DOE_forest Mdl_pr_forest DOE_pr_forest prop_pr_forest...
    prop_pr_DOE   err_flag output t count INPUT RUN result FUNoutput input

% results_satisfy(1:56)=0;
% results_valid(1:56)=0;
%     MissionName='GEP HB';
% MissionName='Lisam210';
% MissionName='3DR Iris';
% MissionName='DJI MATRICE 100';
MissionName='DJI SW 1000';
% MissionName='DJI MG1';
% MissionName='DJI SW 1000 NrNcdp';
%% read missions, constraints and constants
% read mission specifications from a given prescenario
MissionSpecs=xlsread('MissionSpecs.xlsx',MissionName);
Ph=MissionSpecs(1,:);   % the ID of Phases
FT=MissionSpecs(3,:);   % the flight time of each phase
for i=1:length(Ph)
    Dp(i,1)=MissionSpecs(4,i);    % the x (horizontal) displacement
    Dp(i,2)=MissionSpecs(5,i);    % the y (vertical) displacement
    v(i,1)=Dp(i,1)/FT(i);     % the x average velocity, i.e. the cruise speed
    v(i,2)=Dp(i,2)/FT(i);     % the y average velocity
end
m_pl=MissionSpecs(6,:); % the payload mass in each phase
P_pl=MissionSpecs(7,:); % the payload power in each phase
% read constraints from a given scenario
Cstr=xlsread('Constraints.xlsx',MissionName);
Cstr=Cstr(:,1);

Constants=xlsread('Constants.xlsx',MissionName);
C=Constants(1:39,1);              % constants
n_M=length(Constants(42,:));    % number of materials
MP=Constants(42:49, 1:n_M);   % material properties table
MPr=Constants(53:60,1:n_M);   % material bill table (price based on shape)
%% read design parameter options
DesignParameters=xlsread('DesignParameters.xlsx',MissionName);
%read the swithes of design parameters
t=1;
for i=1:6
    if DesignParameters(i,3)==1
        N_r(t)=DesignParameters(i,4);   %read the available number of rotors
        cx(t)=DesignParameters(i,5);    %read the coaxial switch
        N_arm(t)=N_r(t)/(cx(t)+1);      %work out number of arms
        t=t+1;
    end
end
%read number of Li-Po battery cells. The spreadsheet defined the middle
%recommended vale and a range, the below code convert these two numbers
%into an array.
N_c(1)=DesignParameters(8,2)-DesignParameters(9,2);
for i=1:(DesignParameters(9,2)*2)
    N_c(i+1)=N_c(i)+1;
end
for i=1:(DesignParameters(9,2)*2+1)
    if N_c(i)<0
        N_c(i)=0;
    end
end
N_c=nonzeros(N_c)';
%read propeller diameters
t=1;
if DesignParameters(13,2)==0
    d_p=DesignParameters(11,2);
else
    for i=DesignParameters(11,2):DesignParameters(13,2):DesignParameters(12,2)
        d_p(t)=i;
        t=t+1;
    end
end
%read number of propeller blades
t=1;
for i=1:8
    if DesignParameters(i+18,3)==1
        Shape(t)=DesignParameters(i+18,1);
        t=t+1;
    end
end
t=1;
for i=1:5
    if DesignParameters(i+27,3)==1
        Mat(t)=DesignParameters(i+27,1);
        t=t+1;
    end
end
%read body material selection
t=1;
for i=1:5
    if DesignParameters(i+33,3)==1
        Mat_body(t)=DesignParameters(i+33,1);
        t=t+1;
    end
end
%% load regression forest models for motor mass and prop price estimation
load("Motor_forest.mat");
load("DOE_forest.mat");
load("Mdl_pr_forest.mat");
load("DOE_pr_forest.mat");
load("prop_pr_forest.mat");
load("prop_pr_DOE.mat");
%% pre-check of the design parameters, call level 2 process function and constraint check
t=1;
% Cstr(15)=10;
Cstr(17)=200;
Cstr(16)=14;
count=1;
RUN=1;
DsPm(1)=8;       %N_r(i);
DsPm(2)=6;       %N_c(j);
DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
DsPm(4)=15*0.0254;       %d_pm;
DsPm(5)=1;       %Mat_body(l);
DsPm(6)=8;       % N_arm

fun_dl=@level2_NrNcdp;    % set the optimisation function
% Disk Loading is the only variable to be changed
ub_dl=Cstr(17);    % set upper boundary
lb_dl=Cstr(16);    % set lower boundary
DL0=lb_dl+0.001;  % set the starting point
options=optimset('Display','iter','TolFun',1e-6,'TolX',1e-5);
 %% changing DL
 DsPm(4)=19*0.0254;
 x=100;
for i=1:x
    DL(i,t)=(Cstr(17)-Cstr(16))/(x-1)*(i-1)+Cstr(16);
    err_per(i,t)=level2_NrNcdp(DL(i,t));
    output1(i,:)=FUNoutput;
end
t=t+1;
DsPm(4)=13*0.0254;
for i=1:x
    DL(i,t)=(Cstr(17)-Cstr(16))/(x-1)*(i-1)+Cstr(16);
    err_per(i,t)=level2_NrNcdp(DL(i,t));
    output2(i,:)=FUNoutput;
end
t=t+1;
DsPm(4)=23*0.0254;
for i=1:x
    DL(i,t)=(Cstr(17)-Cstr(16))/(x-1)*(i-1)+Cstr(16);
    err_per(i,t)=level2_NrNcdp(DL(i,t));
    output3(i,:)=FUNoutput;
end
t=t+1;
DsPm(4)=19*0.0254;
Cstr(1)=15;     % max speed 15m/s
for i=1:x
    DL(i,t)=(Cstr(17)-Cstr(16))/(x-1)*(i-1)+Cstr(16);
    err_per(i,t)=level2_NrNcdp(DL(i,t));
    output4(i,:)=FUNoutput;
end
t=t+1;
Cstr(1)=25;     % max speed 25m/s
x=30;
Cstr(17)=100;
Cstr(16)=40;
for i=1:x
    DL(i,t)=(Cstr(17)-Cstr(16))/(x-1)*(i-1)+Cstr(16);
    err_per(i,t)=level2_NrNcdp(DL(i,t));
    output5(i,:)=FUNoutput;
end
t=t+1;

Cstr(17)=85;
Cstr(16)=75;

for i=1:x
    DL(i,t)=(Cstr(17)-Cstr(16))/(x-1)*(i-1)+Cstr(16);
    err_per(i,t)=level2_NrNcdp(DL(i,t));
    output6(i,:)=FUNoutput;
end
t=t+1;

%% max speed
for i=1:10
    DL=64+0.2*i;
    err(i)=level2_NrNcdp(DL)
end


%% changing dp and optimise MTOM
for n=7:2:25
    DsPm(4)=n*0.0254;       %d_pm;
    time1=tic;
    [opt_result,Min_err_per,existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);


    
    runtime1=toc(time1);  
    
    output(count,:)=FUNoutput;
    count=count+1;  
    time2=tic;
    [opt_result(count),Min_err_per(count),existflag]=fmincon(fun_dl,DL0,[],[],[],[],lb_dl,ub_dl,[], options);
    runtime2=toc(time2);  
    
    output(count,:)=FUNoutput;
    count=count+1;
%     input=[8,6,0.4826,1,DL];
  fun=@level2_opt_NrNcdp;  
options_l1=optimoptions('ga','PlotFcn',{@gaplotbestf, @gaplotstopping,@gaplotbestindiv},'PopulationSize',10,...
    'MaxGenerations',10,'MaxStallGenerations', 4,'FunctionTolerance',1e-2);
% lb_l1=[Cstr(16),Cstr(14)];    % N_r=x1*2, N_c=x2*2,d_p=13+2*(x3-1), DL, TrW
% ub_l1=[Cstr(17),Cstr(15)];   % lower boundary of variables
%     x0=(lb_l1+ub_l1)/2;
 nonlcon_dl=@level2_NrNcdp_con;
% nonlcon=nonlcon_l1;

    time3=tic;
        [DL0,fval(count),exitflag] =...
            ga(fun,1,[],[],[],[],lb_dl,ub_dl,nonlcon_dl,[],options_l1);
    runtime3=toc(time3);  
        [opt_result(count),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
    runtime4=toc(time3);
    output(count,:)=FUNoutput;
    count=count+1;
end
DsPm(4)=15*0.0254;       %d_pm;
count=count+1;
%% changing Nc and optimise MTOM
for n=2:12
    DsPm(2)=n;       %N_c(j);    
    [opt_result(count),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
    output(count,:)=FUNoutput;
    count=count+1;    
end
DsPm(2)=6;       %N_c(j);
count=count+1;
%% changing Mm/MTOM ratio and optimise MTOM
for n=1:11
     C(30)=0.1+(n-1)*0.02;
    [opt_result(count),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
    output(count,:)=FUNoutput;
    count=count+1;    
end
C(30)=0.1756;      
count=count+1;

%% changing prop FM and optimise MTOM
for n=1:11
    C(6)=0.6+(n-1)*0.02;
    [opt_result(count),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
    output(count,:)=FUNoutput;
    count=count+1;    
end
C(6)=0.68;       %N_c(j);
count=count+1;

%% changing motor lowest efficiency and optimise MTOM
for n=1:11
    C(14)=0.5+(n-1)*0.02;
    [opt_result(count),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
    output(count,:)=FUNoutput;
    count=count+1;    
end
C(6)=0.6;       %N_c(j);
count=count+1;

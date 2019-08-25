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
MissionName='diffconstraints';
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
m_pl=MissionSpecs(6,:); % the payload mass in each phase, uint:g
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
% Cstr(17)=200;
% Cstr(16)=14;
count=1;
RUN=1;
DsPm(1)=N_r;
DsPm(2)=N_c;
DsPm(3)=2;     %N_bl(k);number of blade not considered in this version
DsPm(4)=d_p*0.0254;       %d_pm;
DsPm(5)=1;       %Mat_body(l);
DsPm(6)=N_arm;

fun_dl=@level2_NrNcdp;    % set the optimisation function
% Disk Loading is the only variable to be changed
ub_dl=Cstr(17);    % set upper boundary
lb_dl=Cstr(16);    % set lower boundary
DL0=lb_dl+0.001;  % set the starting point
options=optimset('Display','iter','TolFun',1e-6,'TolX',1e-5);
 %% changing PL
 x=11;
for i=1:x
    m_pl=0.3*(i-1);
    [opt_result(count,:),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
results(count,:)=FUNoutput;
count=count+1;

end
count=count+10;
m_pl=0.5;

 %% changing max speed
 x=9;
for i=1:x
    Cstr(1)=6+3*(i-1);
    [opt_result(count,:),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
results(count,:)=FUNoutput;
count=count+1;
end
count=count+10;
Cstr(1)=15;

 %% changing max endurance
 x=8;
for i=1:x
    FT=(2+4*(i-1))*60;
    [opt_result(count,:),Min_err_per(count),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
results(count,:)=FUNoutput;
count=count+1;
end
count=count+10;
FT=600;

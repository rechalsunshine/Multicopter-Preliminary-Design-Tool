clear all
clc
close all
global Ph FT Dp v m_pl P_pl Cstr C MP cx Shape Mat DesignParameters DsPm MPr ...
    Motor_forest DOE_forest Mdl_pr_forest DOE_pr_forest prop_pr_forest...
    prop_pr_DOE   err_flag output t count INPUT RUN result result_style 

% results_satisfy(1:56)=0;
% results_valid(1:56)=0;
%     MissionName='GEP HB';

MissionName='style';
% MissionName='3DR Iris';
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
Cstr(15)=7;
count=1;
RUN=1;
INPUT(1)=3;             %N_r
INPUT(2)=4;             %N_c
INPUT(3)=1;             %Mat_body
INPUT(4)=10;            %dp
DsPm(1)=INPUT(1)*2;       %N_r(i);
DsPm(2)=INPUT(2);       %N_c(j);
DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
DsPm(4)=INPUT(4)*0.0254;       %d_pm;
DsPm(5)=INPUT(3);       %Mat_body(l);
DsPm(6)=INPUT(1)*2;       % N_arm
fun=@level2_opt_NrNcdp;
options_l1=optimoptions('ga','PlotFcn',{@gaplotbestf, @gaplotstopping,@gaplotbestindiv},'PopulationSize',5,...
    'MaxGenerations',3,'MaxStallGenerations', 2,'FunctionTolerance',1e-5);
lb_l1=[Cstr(16),Cstr(14)];    % N_r=x1*2, N_c=x2*2,d_p=13+2*(x3-1), DL, TrW
ub_l1=[Cstr(17),Cstr(15)];   % lower boundary of variables
%     x0=(lb_l1+ub_l1)/2;
nonlcon_l1=@level2con_NrNcdp;
nonlcon=nonlcon_l1;
%     options=optimset('Display','iter','TolFun',1e-3,'TolX',1e-3,'MaxFunEvals',2000);
%
fun_dl=@level2;    % set the optimisation function
% Disk Loading is the only variable to be changed
ub_dl=Cstr(17);    % set upper boundary
lb_dl=Cstr(16);    % set lower boundary
DL0=(ub_dl+lb_dl)/2;  % set the starting point
options=optimset('Display','iter','TolFun',1e-4,'TolX',1e-4);
    
    dltrw=[0,0];
    min_err(count)=10;
    for i=1:13
        TrW=(Cstr(15)-Cstr(14))/10*(i-1)+Cstr(14);
        % define the design parameters to be used in sub-level functions
        DsPm(7)=TrW;       % TrW
        
        [opt_result,Min_err_per,existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
        if Min_err_per<0.05
            dltrw=[opt_result,TrW];
            break
        elseif Min_err_per<min_err(count)
                min_err(count)=Min_err_per;
                dltrw=[opt_result,TrW];
        end
    end
    %     [depm_opt(count,:),fval(count),exitflag] =...
    %         ga(fun,2,[],[],[],[],lb_l1,ub_l1,nonlcon_l1,[],options_l1);
    if dltrw~=[0,0]
         x0=[INPUT,dltrw];
        result_style(count,:)=level2_opt(x0);
        if result_style(count,19)<0.05
        [depm_opt(count,:),fval(count),exitflag]=fminsearchcon(fun,dltrw,lb_l1,ub_l1,[],[],nonlcon, options);
        dltrw_opt=depm_opt(count,:);
        input1=[INPUT,dltrw_opt];
        result_style(count,:)=level2_opt(input1);
        end
        count=count+1;
    else        
        result_style(count,2:6)=[INPUT(1:2),INPUT(4),INPUT(3)*25.4];
        count=count+1;
    end



Dp_tot=500; % total horizontal displacement 500 m
Dp(1:2,2)=[50,0];
v(1:2,1)=10;      % constant forward speed 10 m/s
funst=@level2_style;
count=2;
options=optimset('Display','iter','TolFun',1e-7,'TolX',1e-7);    

for n=0:5
    Dp(1,1)=100*n;
    Dp(2,1)=Dp_tot-100*n;
    if n==0
        FT(:,1)=[10;50]; 
        v(1,1)=0;
    else
    FT=Dp(:,1)./v(:,1);
    v(1:2,2)=Dp(:,2)./FT;
    end
    if n==5
        v(2,2)=0;
    end
    [opt_result,Min_err_per,existflag]=fminsearchcon(funst,dltrw_opt(1),lb_dl,ub_dl,[],[],[], options);
        As=pi*(DsPm(4)/2)^2;            % disk area of a single rotor
        A=As*DsPm(1);                   % total disk area
        GTOW=A*opt_result;                      % initial GTOW, disk loading times
        % total disk area, unit N. not including payload weight
        GTOM0=GTOW/C(3,1);                % convert weigt to mass, N to kg
        TrW=result_style(1,8)*result_style(1,17)/GTOM0/1000;  
         input2=[INPUT,opt_result,TrW];
        result_style(count,:)=level2_opt(input2);
        count=count+1; 
end

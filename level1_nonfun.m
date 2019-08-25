
clear all
clc
close all
global Ph FT Dp v m_pl P_pl Cstr C MP cx Shape Mat DesignParameters DsPm MPr ...
    Motor_forest DOE_forest Mdl_pr_forest DOE_pr_forest prop_pr_forest...
    prop_pr_DOE   err_flag output err_percent result t count FUNoutput RUN

MissionName='DJI SW 1000';
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
count=1;
RUN=1;
fun_dl=@level2_NrNcdp;    % set the optimisation function
% Disk Loading is the only variable to be changed
ub_dl=Cstr(17);    % set upper boundary
lb_dl=Cstr(16);    % set lower boundary
DL0=lb_dl+0.001;  % set the starting point
options=optimset('Display','iter','TolFun',1e-6,'TolX',1e-5);
    timeoa=tic;
for i=1:length(N_arm)           %loop all the number of rotors
    for j=1:length(N_c)       %loop all the number of battery cells
        %for k=1:length(N_bl)
        for l=1:length(Mat_body)
            for m=1:length(d_p)   %loop all the propeller sizez
                d_pm=d_p(m)*25.4/1000;          % convert propeller diameter from inch to m
DsPm(1)=N_r(i); 
DsPm(2)=N_c(j);
DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
DsPm(4)=d_pm;
DsPm(5)=Mat_body(l);
DsPm(6)=N_arm(i);

    [opt_result(t,:),Min_err_per(t),existflag]=fminsearchcon(fun_dl,DL0,lb_dl,ub_dl,[],[],[], options);
results(t,:)=FUNoutput;
t=t+1;
            end
        end
    end
end

    runtimeoa=toc(timeoa);



% 
% 
% t=1;
% count=1;
% opt_target=0;
% for n=1:length(opt_target)
%     Cstr(18)=opt_target(n) ;
%     
%     % set the optimisation function
%     switch opt_target(n)
%         case 0
%             fun=@level2_opt_MTOM;
%         case 1
%             fun=@level2_opt_Price;
%         case 2
%             fun=@level2_opt_Energy;
%         case 3
%             fun=@level2_opt_EnHo;
%         case 4
%             fun=@level2_opt_EnCru;
%         case 5
%             fun=@level2_opt_Vfw;
%         otherwise
%             fun=@level2_opt_Vas;
%     end
%     
%     options_l1=optimoptions('ga','PlotFcn',{@gaplotbestf, @gaplotstopping,@gaplotbestindiv},'PopulationSize',7,...
%         'MaxGenerations',200,'MaxStallGenerations', 15,'FunctionTolerance',1e-6);
%     
%     intcon=[1,2,3,4];                       % Integer variables
%     lb_l1=[min(N_r),min(N_c),min(Mat_body),min(d_p)];    % N_r=x1*2, N_c=x2*2,d_p=13+2*(x3-1), DL, TrW
%     ub_l1=[max(N_r),max(N_c), max(Mat_body),max(d_p)];   % lower boundary of variables
%     % Disk Loading and thrust weight ratio are the two variables
%     %                     ub_l1=[Cstr(17),Cstr(15)];    % set upper boundary
%     %                     lb_l1=[Cstr(16),Cstr(14)];    % set lower boundary
%     % %                     input0=[(Cstr(17)+Cstr(16))/2,Cstr(15)*0.2];  % set the starting point
% %     nonlcon_l1=@level2con;
%     timeoa=tic;
%     [depm_opt(n,:),fval,exitflag] =...
%         ga(fun,4,[],[],[],[],lb_l1,ub_l1,[],intcon,options_l1);
%     % [depm_opt(t,:),fval,exitflag] =...
%     %     ga(fun,2,[],[],[],[],lb_l1,ub_l1,nonlcon_l1,[],options_l1);
%     runtimeoa=toc(timeoa);
%     input1=depm_opt(n,:);
% %     % define the design parameters to be used in sub-level functions
% %                     DsPm(1)=depm_opt(1)*2;       %N_r(i);
% %                     DsPm(2)=depm_opt(2);       %N_c(j);
% %                     DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
% %                     DsPm(4)=depm_opt(4)*0.0254;       %d_pm;
% %                     DsPm(5)=depm_opt(3);       %Mat_body(l);
% %                     DsPm(6)=depm_opt(1)*2;       % N_arm
% %                     DsPm(7)=depm_opt(6);       % TrW
% %                     DL1=depm_opt(5);
%     performance_opt(n,:)=level2_opt_MTOM(input1);
%     results_opt(n,:)=FUNoutput;
% end
%     result_valid=result;
%     
%     n=2;
%     fun=@level2_opt_MTOMwrap;
%   options_l1=optimoptions('ga','PlotFcn',{@gaplotbestf, @gaplotstopping,@gaplotbestindiv},'PopulationSize',30,...
%         'MaxGenerations',200,'MaxStallGenerations', 15,'FunctionTolerance',1e-6);
%     
%     intcon=[1,2,3,4];                       % Integer variables
%     lb_l1=[min(N_r),min(N_c),min(Mat_body),min(d_p),Cstr(16)];    % N_r=x1*2, N_c=x2*2,d_p=13+2*(x3-1), DL, TrW
%     ub_l1=[max(N_r),max(N_c), max(Mat_body),max(d_p),Cstr(17)];   % lower boundary of variables
%     % Disk Loading and thrust weight ratio are the two variables
%     %                     ub_l1=[Cstr(17),Cstr(15)];    % set upper boundary
%     %                     lb_l1=[Cstr(16),Cstr(14)];    % set lower boundary
%     % %                     input0=[(Cstr(17)+Cstr(16))/2,Cstr(15)*0.2];  % set the starting point
%     nonlcon_l1=@level2con;
%     timeoa1=tic;
%     [input2,fval,exitflag] =...
%         ga(fun,5,[],[],[],[],lb_l1,ub_l1,nonlcon_l1,intcon,options_l1);
%     % [depm_opt(t,:),fval,exitflag] =...
%     %     ga(fun,2,[],[],[],[],lb_l1,ub_l1,nonlcon_l1,[],options_l1);
%     runtimeoa1=toc(timeoa1);
%         performance_opt(n,:)=level2_opt_MTOMwrap(input2);
%     results_opt(n,:)=FUNoutput;

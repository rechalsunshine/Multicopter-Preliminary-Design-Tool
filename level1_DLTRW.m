clear all
clc
close all
global Ph FT Dp v m_pl P_pl Cstr C MP cx Shape Mat DesignParameters DsPm MPr ...
    Motor_forest DOE_forest Mdl_pr_forest DOE_pr_forest prop_pr_forest...
    prop_pr_DOE   err_flag output

results_satisfy(1:56)=0;
results_valid(1:56)=0;
%     MissionName='GEP HB';

% MissionName='DJI SW 1000';
MissionName='3DR Iris';
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
i=1;
j=1;
% for a=1:length(N_arm)           %loop all the number of rotors
%     for b=1:length(N_c)       %loop all the number of battery cells
%         %for k=1:length(N_bl)
%         for c=1:length(Mat_body)
%             for d=1:length(d_p)   %loop all the propeller sizez
%                 d_pm=d_p(d)*25.4/1000;          % convert propeller diameter from inch to m
n=20;
for i=1:n
    for j=1:n
        DL(i)=(Cstr(17)-Cstr(16))/(n-1)*(i-1)+Cstr(16);
        TrW(j)=(Cstr(15)-Cstr(14))/(n-1)*(j-1)+Cstr(14);
        % scan the range of thrust weight ratio, divide it into 10 attempts.
        
        % define the design parameters to be used in sub-level functions
        DsPm(1)=4;          %N_r(a);
        DsPm(2)=3;          %N_c(b);
        DsPm(3)=2;     %N_bl(k);number of blade not considered in this version
        DsPm(4)=0.2667;     %d_pm;        %d_pm;
        DsPm(5)=2;          %Mat_body(c);
        DsPm(6)=4;          %N_r(a);      % N_arm(i);
        DsPm(7)=TrW(j);
%         DL=55.8376;     
        
        err_percent(i,j)=level2(DL(i));
                
        results_all_dltrw(t,:)=output;    % record all the valid
        % results, but some of them may not satisfy the
        % mission requirements.
        t=t+1;
        
%         if sum(err_flag)==0
%             % if the results pass the constraints check of
%             % the mission requiremnts, then record the
%             % result and parameters.
%             results_satisfy_dltrw(count,:)=output;
%             
%         end
        fprintf('loop %d, %d. error percentage is: %f. \n', i,j,err_percent(i,j));
%     end
% end
count=count+1;
%             end
%         end        
    end
end

% %%%%%  optimisation process %%%%
% %1. GA for rough optimisation
% optionsga=optimoptions('ga','PopulationSize',20,'MaxGenerations',40,...
%     'MaxStallGenerations', 10,'FunctionTolerance',1e-5,'Display','off');
% intcon=[1,2,3];                       % Integer variables
% nonlconga=@level3_armcons;
% lbga=[2,2, 1,Cstr(16),Cstr(14)];   % lower boundary of variables
% ubga=[4,4, 3,Cstr(17),Cstr(15)];    % N_r=x1*2, N_c=x2*2,d_p=13+2*(x3-1), DL, TrW
% 
% [depm_opt(1:5),MTOM_min,exitflag] =...
%     ga(@level2_opt,5,[],[],[],[],lbga,ubga,nonlconga,intcon,optionsga);
% create a plot of the correlations
% figure;
% surf(TrW,DL,err_percent)
% legend('MTOM error');
% ylabel('Disk Loading'); xlabel('Thrust to weight ratio');

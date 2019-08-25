clear all
clc
close all
global Constants MissionSpecs MaterialProperties Shape Material...
    f_sep cdl SF PlatformParams output s Ti
%MTOM=1.4;
%MissionName='MGatti_Hover';
MissionName='Default';

%read givne platform criteria

PlatformParams=xlsread('database.xlsm',MissionName);

%read mission parameters from a given scenario
MissionSpecs=xlsread('Missions.xlsx',MissionName);

%read fixed parameters from a given scenario
Constants=xlsread('Constants.xlsx','higherefficiency');

%read material database
for i=1:5
    MaterialProperties(1:7,i)=xlsread('Materials.xlsx',i);
end

%read frame preferences
DesignParameters=xlsread('DesignParameters.xlsx',MissionName);
f_sep=DesignParameters(19,3);              % ??? allowed minimum frequency separation
cdl=DesignParameters(20,3);                  % ??? the coefficient of deflection limit
SF=DesignParameters(21,3);                   % ??? yield strength safety factor

t=[1,2,3,4,5,16,28];
r=[5,8,11,14,17,20,23];
s=1;
c=1;
for i=1:length(t)   % loop to go through all the platforms
%     i=2;
    Ti=1;
    % read the given platform parameters
    for j=1:3
        MissionSpecs(3,j)=PlatformParams(t(i),8)/1000;     %payload, unit: kg
    end
    MissionSpecs(7,1)=PlatformParams(t(i),13);     %max ascent velocity, unit: m/s
    MissionSpecs(7,3)=PlatformParams(t(i),14);     %max descent velocity, unit: m/s
    MissionSpecs(2,1)=MissionSpecs(5,1)/MissionSpecs(7,1);   % take-off phase flight time, unit: s
    MissionSpecs(2,3)=MissionSpecs(5,3)/MissionSpecs(7,3);   % take-off phase flight time, unit: s
    MissionSpecs(2,2)=PlatformParams(t(i),11)*60;     % hovering time, unit: s
    MissionSpecs(6,2)=PlatformParams(t(i),12);     % maximum horizontal velocity, unit: m/s
    
    Shape=PlatformParams(t(i),27);     % shape code
    Material=PlatformParams(t(i),20);     % material code
    
    Constants(1,1)=PlatformParams(t(i),38);     % Disk Loading, unit: N/m2
    Constants(2,1)=PlatformParams(t(i),34);     % number of blade
    Constants(4,1)=PlatformParams(t(i),15);     % number of rotors
    Constants(8,1)=PlatformParams(t(i),46);     % number of battery cells
    Constants(19,1)=PlatformParams(t(i),21);     % Drag area, unit: m2
    Constants(21,1)=PlatformParams(t(i),5)/1000;     % initial MTOM, unit: kg
    Constants(23,1)=PlatformParams(t(i),40);     % propeller gap/diameter ratio
    Constants(24,1)=PlatformParams(t(i),22);     % is coaxial? bool
    
    output(6+s)=PlatformParams(t(i),11);     % hovering endurence, unit: min
    output(7+s)=PlatformParams(t(i),12);     % max horizontal speed
    output(8+s)=PlatformParams(t(i),13);     %max ascent velocity, unit: m/s
    output(9+s)=PlatformParams(t(i),14);     %max descent velocity, unit: m/s
    output(10+s)=PlatformParams(t(i),15);     % number of rotors
    output(16+s)=PlatformParams(t(i),21);     % Drag area, unit: m2
    output(32+s)=Constants(1,1);     % Disk Loading, unit: N/m2
    output(34+s)=Constants(23,1);     % propeller gap/diameter ratio
    output(52+s)=Constants(8,1);     % number of battery cells

    
    
%     options=optimoptions('fmincon','Algorithm','sqp-legacy','Display','iter-detailed','CheckGradients', true,'OptimalityTolerance',...
%         1e-8,'StepTolerance', 1e-8,'MaxFunctionEvaluations',...
%         3000, 'MaxIterations',1000); % options for fmincon
    

options=optimset('Display','iter','TolFun',1e-5,'TolX',1e-5);
    MTOM0=Constants(21,1);         % read the initial MTOM, actually the real MTOM of the product
    %MTOM0=3;
    ub=1.5*MTOM0;
    lb=0.1*MTOM0;
   % nonlconf=@MTOMConf;
    timeoa=tic;
%     [est_MTOM1,Min_err1]=fminsearchbnd(@MTOMworkoutFrameGA,MTOM0,lb,ub,options1);
%     [est_MTOM2,Min_err2]=fminsearchbnd(@MTOMworkoutFrameGA,est_MTOM1,lb,ub,options1);
    [est_MTOM,Min_err,existflag]=fminsearchbnd(@MTOMworkoutFrameGA_fix,MTOM0,lb,ub,options);
    runtimeoa(c,1)=toc(timeoa);
    c=c+1;
%     xlrangert=['D',num2str(t+4), ':D',num2str(t+4)];
%     xlswrite('SimulationTrack.xlsx',runtimeoa(t,1),MissionName,xlrangert);
%     
%     output(2)=est_MTOM*1000; % optimised MTOM, unit: g
%     output(3)=Min_err*1000; %  MTOM error, unit: g
   
    xlrange=['C',num2str(r(i)), ':BM',num2str(r(i))];
    xlswrite('SimulationTrack.xlsx',output,'BOMoptimisation_m&b',xlrange);
end
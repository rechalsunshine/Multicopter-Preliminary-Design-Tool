clear all
clc
close all
global Ph FT Dp v m_pl P_pl Cstr C MP cx Shape Mat DesignParameters DsPm MPr ...
    Motor_forest DOE_forest Mdl_pr_forest DOE_pr_forest prop_pr_forest...
    prop_pr_DOE   err_flag output t count INPUT RUN result FUNoutput MidPara Load  missionerr 

%     MissionName='GEP HB';
% MissionName='Lisam210';
MissionName='3DR Iris';
% MissionName='DJI MATRICE 100';
% MissionName='DJI SW 1000';
% MissionName='DJI MG1';
%% read missions, constraints and constants
% read mission specifications from a given prescenario
% MissionSpecs=xlsread('MissionSpecs.xlsx',MissionName);
% Ph=MissionSpecs(1,:);   % the ID of Phases
% FT=MissionSpecs(3,:);   % the flight time of each phase
% for i=1:length(Ph)
%     Dp(i,1)=MissionSpecs(4,i);    % the x (horizontal) displacement
%     Dp(i,2)=MissionSpecs(5,i);    % the y (vertical) displacement
%     v(i,1)=Dp(i,1)/FT(i);     % the x average velocity, i.e. the cruise speed
%     v(i,2)=Dp(i,2)/FT(i);     % the y average velocity
% end
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
% DsPm(1)=INPUT(1)*2;       %N_r(i);
% DsPm(2)=INPUT(2);       %N_c(j);
% DsPm(3)=0;     %N_bl(k);number of blade not considered in this version
% DsPm(4)=INPUT(4)*0.0254;       %d_pm;
% DsPm(5)=INPUT(3);       %Mat_body(l);
% DsPm(6)=INPUT(1)*2;       % N_arm

DsPm=[4,4,0,7*0.0254,1,4];
DL=68.52518;

err_flag(1:9)=0;     % clear the cache
missionerr=[0,0,0];
err_percent=0;
FUNoutput(1:63)=1;

%% fore calculation
% for notion: DsPm refers to design parameters: (1)N_r,number of rotors,
% (2) N_c, number of battery cells, (3) N_bl, number of propeller blades,
% (4) d_p, propeller diameter, unit,m, (5) body material flag, (6) N_arm,
% number of arms, (7) TrW, thrust weight ratio
As=pi*(DsPm(4)/2)^2;            % disk area of a single rotor
A=As*DsPm(1);                   % total disk area
GTOW=A*DL;                      % initial GTOW, disk loading times
% total disk area, unit N. not including payload weight
GTOM0=GTOW/C(3);                % convert weigt to mass, N to kg
MTOM=GTOM0+max(m_pl);   % maximum take off mass, including max payload mass
MTOW=MTOM*C(3);     % max take off weight, including max payload
U=C(1)*DsPm(2);     % total voltage equals cell voltage times number of cells
m_av=Cstr(13);      % avionics mass and power
P_av=Cstr(12);
% central body and landing gear
L_a=DsPm(4)*(1+C(23))/(2*sin(pi/DsPm(6)));  % the arm length of given propeller
% it is defined as from the centre of body to the centre of a motor. By
% defining the gap ratio between two adjacent propellers, given the
% propeller diameter, the arm length can be determined.
L_d=2*L_a;        % diagonal distance is considered twice of the arm length
% define some middle parameters to be used by level3 functions
MidPara(:,1)=[As;0;GTOW;0;MTOM];
% MidParameters={As,'As-Disk Area of one rotor';A1,'A1-topview area';GTOW,'GTOW-gross takeoff weight';...
%     v_h,'v_h-hovering induced velocity';MTOM,'max takeoff mass'};
%% propeller estimation
[m_p,pr_p]=level3_prop(DsPm(4));     %% call a level 3 funciton to calculate the mass and price of a single propeller
%% central body estimation
[m_body,pr_body,A_c, d_body, th_body]=level3_body(L_d);  % mass, price, and the central frame area
%% Landing Gear Constraints
H_lg_max=0.2829*L_d;                    % maximum allowed landing gear height
m_lg=Cstr(21);                          % Landing gear mass
%% Drag calculation from max straight forward condition
% the calculated top area should be within a reseanable range, here
% consider 90% of the total disk area+central part area (10% overlap).
A1=(A+A_c)*0.9;
A2=A1/C(27,1);  % the front area, it is normally the body+payload+landing gear
% projection area, it is much smaller than the top area, according to
% database conclusion, there is an average number
MidPara(2,1)=A1;
%% Hovering condition--with maximum payload
var_h=[0,0,MTOM];
[Ts_h,n_hp,Psp_hp]=level3_general(var_h);
Psp_mot=Psp_hp*Cstr(14);      % when considering maximum thrust mode,
% also need to multiply by the minimum thrust to weight ratio to cover manueverbility
Ts_max(1)=Ts_h*Cstr(14);
Ps_max(1)=Psp_mot;
n_max(1)=n_hp;      % record candidate one for maximum motor output

%%
v_const=10;
step=100;
X=1000;
Y=3000;
n=X/step+1;
T_tot=(X+Y)/v_const;

for i=1:n
    FT(i,2)=(X-(i-1)*step)/v_const;
    FT(i,1)=T_tot-FT(i,2);
    v(i,:)=[(i-1)*step/FT(i,1),Y/FT(i,1),v_const,0];

    %phase1:
    var1=[v(i,1:2),MTOM];
    [Ts_max(i,1),n_f_max(i,1),Ps(i,1)]=level3_general(var1);    
        Psp(i,1)=Ps(i,1)*C(24)/ C(7) ;
        
        %phase2:
    var2=[v(i,3:4),MTOM];
    [Ts_max(i,2),n_f_max(i,2),Ps(i,2)]=level3_general(var2);    
        Psp(i,2)=Ps(i,2)*C(24)/ C(7) ;
 
        E(i,1:2)=Psp(i,:).*FT(i,:)./3600;
        E(i,3)=E(i,1)+E(i,2);
  
end
out=[FT,v(:,1:2),Psp,E];
for i=12:21
        FT(i,2)=(i-11)*150/10;
    FT(i,1)=300-FT(i,2);
    v(:,:,i)=[1500/FT(i,1),(1500-(i-11)*150)/FT(i,1);0,10];

    %phase1:
    var1=[v(1,:,i),MTOM];
    [Ts_max(i,1),n_f_max(i,1),Ps(i,1)]=level3_general(var1);    
        Psp(i,1)=Ps(i,1)*C(24)/ C(7) ;
        
        %phase2:
    var2=[v(2,:,i),MTOM];
    [Ts_max(i,2),n_f_max(i,2),Ps(i,2)]=level3_general(var2);    
        Psp(i,2)=Ps(i,2)*C(24)/ C(7) ;
 
        E(i,1:2)=Psp(i,:).*FT(i,:)./60;
        E(i,3)=E(i,1)+E(i,2);
end
for i=1:21
    vfa(i,:)=v(1,:,i);
end
for i=1:5
    v(i,1)=1e-5*10^i;
    v(i,2)=10-v(i,1);
    var=[v(i,:),MTOM];
    [Ts_max(i,1),n_f_max(i,1),Ps(i,1)]=level3_general(var); 
end

%% max straight forward flight --with no payload
if Cstr(1)==0
    v_f_max=max(v(:,1));
else
    v_f_max=Cstr(1);     % get maximum forward speed ;
end
var_fmax=[v_f_max,0,GTOM0];
[Ts_max(2),n_f_max,Ps_f_max]=level3_general(var_fmax);
Psp_f_max=Ps_f_max*C(24)/ C(9) ;               % maximum output power happens at extreme
% conditions, the prop efficiency is considered to be the allowed
% worst, 60%. Also considers the average fuselage downwash for multicopters
Ps_max(2)=Psp_f_max;
n_max(2)=n_f_max;      % record candidate two for maximum motor output
% power and rotation speed
%% maximum vertical  ascent condition---no payload
if Cstr(2)==0
    v_a_max=max(v(:,2));
else
    v_a_max=Cstr(2);
end
var_amax=[0,v_a_max,GTOM0];
[Ts_max(3),n_a_max,Ps_a_max]=level3_general(var_amax);
Psp_a_max=Ps_a_max*C(24)/ C(9);               % maximum output power happens at extreme
% conditions, the prop efficiency is considered to be the allowed
% worst, 60%. Also considers the average fuselage downwash for multicopters
Ps_max(3)=Psp_a_max;
n_max(3)=n_a_max;      % record candidate three for maximum motor output
% power and rotation speed
%% extreme condition analysis
% this condition is to work out the power and current requirement uncer the
% extreme performance, it will then estimate the mass and price of motors,
% ESCs and arms
[Ps_in_max, I]=max(Ps_max);   % find the max power and it's condition
n_max_req=n_max(I);   % rotation speed of the rotor, unit: rev/s
Ts_max_req=Ts_max(I); % maximum thrust, unit: N
if Ts_max_req<max(Ts_max)
    err_flag(3)=1;     % error: max power and max thrust are not the same condition
    Ts_max_req=max(Ts_max);
end
TrW=Ts_max_req*DsPm(1)/MTOW;        % expected TrW
if (TrW<Cstr(14)) || (TrW>Cstr(15))
    err_flag(4)=1;         % TrW not in the appropriate range
end
%Mach0.7 tip velocity check
omega_max=n_max_req*2*pi;
v_tip=omega_max*DsPm(4)/2;
if v_tip/C(4)>C(11)
    err_flag(9)=1;  % the prop tip velocity Mach number is exceeded, possible
    % noise and efficiency drop, should avoid.
    %     FUNoutput(19)=10;
    %     FUNoutput(18)=100000;
    return;
end

% kv_min=ceil(n_max_req*60*C(29)/(U*10))*10;      
kv_min=n_max_req*60*C(29)/U;
% kv is rpm per volt, so the rotation speed needs
% a unit conversion. The kv calculated from the minimum rpm is the minimum
% requirement. The actual kv mush be larger than this number, normally, there
% is a ratio between the noted (unloaded) kv vs the full-loaded kv: 1.52

% There are three parameters used to estimate the mass of the motor: 1. minimum kv,
% 2. maximum output power, 3. maximum voltage
param_motor(1)=kv_min;
param_motor(2)=Ps_in_max;
param_motor(3)=U;
[m_mot,pr_mot]=level3_motor(param_motor);
%% ESC build up and mass estimation
P_mt_in_max=Ps_in_max/C(14);    % C(14) is the lowest efficiency of motor
% in the extreme operation condition. !!!!need verify
P_esc_in_max=P_mt_in_max/C(18); % C(18) is the lowest ESC efficiency
I_dc_max=P_esc_in_max/U;           % max current through a single ESC

% There are two parameters used to estimate the mass and price of the
% ESC: 1. maximum input current, 2. the arm length
param_esc(1)=I_dc_max;
param_esc(2)=L_a;
[m_esc,pr_esc,m_esc_chip]=level3_esc(param_esc);
%%  efficient condition
% normally it is hover or forward cruise, but it could also be forward
% and climbing at the same time in a low speed.But no matter what, the
% movemenets in the mission specification are all considered in the
% efficient mode.
for i=1:length(Ph)
    GTOM_ph=GTOM0+m_pl(i);
    var=[v(i,1),v(i,2),GTOM_ph];    % variables for general function
    [~,~,Psp_g]=level3_general(var);
    Ps_g_mot=Psp_g*C(24)/C(7);     % the motor output power in other conditions.
    % consider maximum prop efficiency here.also consider fuselage downwash
    Ps_g=Ps_g_mot/(C(12)*C(17));   % the ESC input power, consider in efficient mode
    % all the components are working at the highest efficiency
    P_ph(i)=(Ps_g*DsPm(1)+P_pl(i)+P_av)/C(20); % total power equals to the
    % power consumed by all the motors, the payload and the avionics. unit: W
    % consider the battery internal loss, apply the efficiency factor
    E_ph(i)=P_ph(i)*FT(i)/3600;     % energy consumed in this phase, unit: W*h
end
E_tot=sum(E_ph);        % total energy
FT_tot=sum(FT);         % total flight time, unit: s
%% calculate the possible endurance
% .1. Hovering with average payload
GTOM_e=GTOM0+Cstr(6);   % here not consider a specific mission, but use an
% average payload mass. Below, uses average payload power.
var_he=[0,0,GTOM_e];
[~,n_he,Psp_he]=level3_general(var_he);
Ps_he_mot=Psp_he*C(24)/C(7);     % the motor output power in other conditions.
% consider maximum prop efficiency here.also consider fuselage downwash
Ps_he=Ps_he_mot/(C(12)*C(17));   % the ESC input power, consider in efficient mode
% all the components are working at the highest efficiency
P_he=(Ps_he*DsPm(1)+Cstr(5)+P_av)/C(20);
FT_he=floor(E_tot/P_he*3600);  % hovering endurance, unit: s
% .2. Hovering without payload
var_he2=[0,0,GTOM0];
[~,n_he2,Psp_he2]=level3_general(var_he2);
Ps_he_mot2=Psp_he2*C(24)/C(7);     % the motor output power in other conditions.
% consider maximum prop efficiency here.also consider fuselage downwash
Ps_he2=Ps_he_mot2/(C(12)*C(17));   % the ESC input power, consider in efficient mode
% all the components are working at the highest efficiency
P_he2=(Ps_he2*DsPm(1)+Cstr(5)+P_av)/C(20);
FT_he2=floor(E_tot/P_he2*3600);  % hovering endurance, unit: s
% .3. Cruising with average payload
var_ce=[Cstr(7),0,GTOM_e];
[~,n_ce,Psp_ce]=level3_general(var_ce);
Ps_ce_mot=Psp_ce*C(24)/C(7);     % the motor output power in other conditions.
% consider maximum prop efficiency here.also consider fuselage downwash
Ps_ce=Ps_ce_mot/(C(12)*C(17));   % the ESC input power, consider in efficient mode
% all the components are working at the highest efficiency
P_ce=(Ps_ce*DsPm(1)+Cstr(5)+P_av)/C(20);
FT_ce=floor(E_tot/P_ce*3600);  % cruising endurance, unit: s
%% Battery calculation and mass estimation
%E_max=max([E_tot,Ea_h,Ea_c]);   % the battery should be enough to provide
% energy in all the conditions (meet the mission requirement as well
% the general endurance constraints)
E_b_req=E_tot/C(21); %  consider battery discharge depth
[m_b,pr_b]=level3_bat(E_b_req) ;    % estimate the battery mass and price from the energy view

Cb=E_b_req/U*1000; %calculate the capacity of the battery, unit: mAh

I_b_max=(P_esc_in_max*DsPm(1)+P_av+max(P_pl))/(C(20)*U);   %work out the maximum curernt in the system to calculate C rate of battery
C_rate=ceil(I_b_max*1000/Cb/5)*5;    %calculate the required C rate (discharge speed)
%% arm optimisation and mass estimation
% La_min=(dt+Dp)/2/sind(Alpha);     % La is the minimum length of arm

Load(1)=Ts_max_req;              % the maximum thrust of a single rotor, unit N
Load(2)=max(n_max);       % the minimum required rotation speed at maximum thrust, unit rev/s
Load(3)=max(Ps_max)/omega_max;     % the maximum torque, unit N.m=W/(rad/s)
Load(4)=m_mot;              % the motormass unit kg
Load(5)=n_he;               % hovering motor rotation speed
Load(6)=n_ce;               % cruising motor rotation speed
if Cstr(18)==1
    [m_arm, pr_arm,opt_geo]=level3_arm_propt(L_a);
else
    [m_arm,pr_arm, opt_geo]=level3_arm_mopt(L_a);
end
%% final add up
 M_tot=DsPm(1)*(m_mot+m_esc+m_p)+DsPm(6)*m_arm+m_b+m_body+m_av+m_lg;
M_err=abs(M_tot-GTOM0);
err_percent=M_err/GTOM0;

Pr_tot=DsPm(1)*(pr_mot+pr_esc+pr_p)+DsPm(6)*pr_arm+pr_b+pr_body;
% price estimation not including avionics and payload. Since they are
% prepared in advance.

%% performance constraints check
if L_d>Cstr(9) && Cstr(9)~=0
    err_size=1; % diagonal distance exceed max size constraint
    err_flag(1)=err_size;
end
if H_lg_max<Cstr(20)
    errC_lg_h=1;     % if this flag is on, it means this solution (if MTOM
    % is valid)is not acceptable due to a "Constraint mismatch of landinggear height"
    % the selected landing gear is too tall for this solution's frame size
    err_flag(2)=errC_lg_h;
end
% if v_f_max<Cstr(1)
%     err_max_vf=1;   % the error flag of max forward speed not meet the mission requirment
%     err_flag(3)=err_max_vf;
% end
% if v_a_max<Cstr(2)
%     err_max_va=1;   % error flag of max ascent speed not meet the mission requirement
%     err_flag(4)=err_max_va;
% end
if FT_he<Cstr(3)
    err_FT_h=1;     % error flag of hovering endurance not meet the mission requirement
    err_flag(5)=err_FT_h;
end
if FT_ce<Cstr(4)
    err_FT_c=1;     % error flag of cruising endurance not meet the mission requirement
    err_flag(6)=err_FT_c;
end
if M_tot+max(m_pl)>Cstr(8) && Cstr(8)~=0
    err_MTOM=1;
    err_flag(7)=err_MTOM;
end
% err_flag(8) is an mistake of vertical descent speed, aircraft will fall
% into vortex ring state.
% err_flag(9) is an mistake that the propeller rotate too fast and the tip
% velocity exceeds the critical Mach number (Mach0.7).
% err_flag(3) max power and max thrust are not the same condition
% err_flag(4) is the indicator of thrust weight ratio out of range
%% mission constraints check
if v_f_max<max(v(:,1))
    missionerr(2)=1;   % the error flag of max forward speed not meet the mission requirment
end
if v_a_max<max(v(:,2))
    missionerr(3)=1;   % error flag of max ascent speed not meet the mission requirement
end
%% record the performance and design parameters of this multicopter

FUNoutput=[Cstr(18),DsPm(1:3),DsPm(4)*1000,DsPm(5:6),TrW,DL,opt_geo(1),...
    opt_geo(2:5).*1000,opt_geo(6),L_d*1000,d_body*1000,th_body*1000,GTOM0*1000,M_tot*1000,...
    err_percent, Pr_tot,E_tot,FT_tot/60,FT_he/60,FT_he2/60,FT_ce/60,v_f_max,...
    v_a_max,m_mot*1000,m_esc*1000,m_esc_chip*1000,m_p*1000,m_arm*1000,m_b*1000,m_body*1000,...
    m_av*1000,m_lg*1000, pr_mot,pr_esc,pr_p,pr_arm,pr_b,pr_body,P_av,Ps_in_max,Ts_max_req...
    Psp_mot,omega_max,kv_min,I_dc_max,Cb,C_rate,H_lg_max*1000,err_flag,missionerr];

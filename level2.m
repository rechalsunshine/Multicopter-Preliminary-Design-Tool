function err_percent=level2(DL)
%timewo=tic;
global Ph FT v Cstr C m_pl P_pl MidPara DsPm output Load err_flag missionerr   
err_flag=[0,0,0,0,0,0,0,0,0];     % clear the cache
output=0;
missionerr=[0,0,0];
err_percent=0;

%% fore calculation
% for notion: DsPm refers to design parameters: (1)N_r,number of rotors,
% (2) N_c, number of battery cells, (3) N_bl, number of propeller blades,
% (4) d_p, propeller diameter, unit,m, (5) body material flag, (6) N_arm,
% number of arms, (7) TrW, thrust weight ratio
As=pi*(DsPm(4)/2)^2;            % disk area of a single rotor
A=As*DsPm(1);                   % total disk area
GTOW=A*DL;                      % initial GTOW, disk loading times
% total disk area, unit N. not including payload weight
GTOM0=GTOW/C(3,1);                % convert weigt to mass, N to kg
MTOM=GTOM0+max(m_pl);   % maximum take off mass, including max payload mass
U=C(1)*DsPm(2);     % total voltage equals cell voltage times number of cells
m_av=Cstr(13);      % avionics mass and power
P_av=Cstr(12);
m_LG=Cstr(21);
% central body and landing gear
L_a=DsPm(4)*(1+C(23))/(2*sin(pi/DsPm(6)));  % the arm length of given propeller
% it is defined as from the centre of body to the centre of a motor. By
% defining the gap ratio between two adjacent propellers, given the
% propeller diameter, the arm length can be determined.
L_d=2*L_a;        % diagonal distance is considered twice of the arm length
TrW=DsPm(7);      % Thrust weight ratio
% define some middle parameters to be used by level3 functions
MidPara(:,1)=[As;0;GTOW;0;MTOM];
% MidParameters={As,'As-Disk Area of one rotor';A1,'A1-topview area';GTOW,'GTOW-gross takeoff weight';...
%     v_h,'v_h-hovering induced velocity';MTOM,'max takeoff mass'};
%% propeller estimation
[m_p,pr_p]=level3_prop(DsPm(4));     %% call a level 3 funciton to calculate the mass and price of a single propeller
%% central body estimation
[m_body,pr_body,A_c]=level3_body(L_d);  % mass, price, and the central frame area
%% Landing Gear Constraints
H_lg_max=0.2829*L_d;                    % maximum allowed landing gear height
%% Drag calculation from max straight forward condition
% the calculated top area should be within a reseanable range, here
% consider 90% of the total disk area+central part area (10% overlap).
A1=(A+A_c)*0.9;
A2=A1/C(27,1);  % the front area, it is normally the body+payload+landing gear
% projection area, it is much smaller than the top area, according to
% database conclusion, there is an average number
MidPara(2,1)=A1;
%% Hovering condition--with maximum payload
T_hp=MTOM*C(3,1);      % equilibrium
Ts_h=T_hp/DsPm(1);               % thrust of a single rotor
n_hp=sqrt(Ts_h/(C(5,1)*C(2,1)*DsPm(4)^4));   % rotation speed of the rotor, unit: rev/s
v_hp=sqrt(Ts_h/(2*C(2,1)*As))  ;    % induced velocity in hover
Ps_hp=v_hp*Ts_h;                  % required hovering power of a single rotor
Psp_hp=Ps_hp*C(24)/C(6);            % prop input power in hover, fixed by FM
% to consider the profile power, and fuselage downwash factor to
% consider the frame area under the propeller.
Psp_mot=Psp_hp*Cstr(14);      % when considering maximum thrust mode,
% also need to multiply by the minimum thrust to weight ratio to cover manueverbility
Ps_max(1)=Psp_mot;
n_max(1)=n_hp;      % record candidate one for maximum motor output
%% max straight forward flight --with no payload
T_max=GTOW*TrW;                 % weight times thrust weight ratio, unit: N
Ts_max=T_max/DsPm(1);           % thrust of a single rotor

beta_f=acos(1/TrW);               % maximum tilt angle,worked out from equilibrium, unit: rad
% beta_f_d=beta_f*180/pi;             % convert rad to degree
alpha_f=beta_f;                     % propeller angle of attack, the angle between
%comming air flow and the propeller disk axis, equals to aircraft tilt
%angle at this straight forward situation. unit: rad
gama_f=0;                         % climb angle, between the incoming wind and horizontal
%direction, equals to vertical speed over forward speed, unit: rad
D=T_max*sin(alpha_f);             % equilibrium
A=A1*sin(alpha_f)+A2*cos(alpha_f);  % the total projected area
q=D/(A*C(26));              % get the dynamic pressure from D=q*Cd*A
v_f_max=sqrt(2*q/C(2));     % get maximum forward speed from q=rho*v_f^2/2;
v_in_f=v_f_max*sin(alpha_f);   % the wind speed that is perpendiculer to the disk plane
var_vi=[alpha_f, v_f_max, Ts_max];
v_i_f=level3_vi(var_vi);        % induced velocity
Ps_f_max=Ts_max*(v_in_f+v_i_f);   % required output power for a single rotor
Psp_f_max=Ps_f_max*C(24)/ C(9) ;               % maximum output power happens at extreme
% conditions, the prop efficiency is considered to be the allowed
% worst, 60%. Also considers the average fuselage downwash for multicopters
n_f=(v_in_f+v_i_f)/(C(8)*C(10)*DsPm(4));   % rotation speed of the rotor, unit: rev/s
% J=v/(n*D), J advance ratio, v incoming flow velocity, D propeller diameter
% when operated at extreme conditions

Ps_max(2)=Psp_f_max;
n_max(2)=n_f;      % record candidate two for maximum motor output
% power and rotation speed
%% maximum vertical  ascent condition---no payload
D_a_max=T_max-GTOW;         % equilibrium when aircraft ascending with max vertical speed
v_a_max=sqrt(D_a_max*2/(C(2)*C(26)*A1));        % max ascending speed
var_vi=[pi/2, v_a_max, Ts_max];
v_i_a=level3_vi(var_vi);        % induced velocity
Ps_a_max=Ts_max*(v_a_max+v_i_a);     % call level3 process to calculate the
%power requirement of a single rotor in ascending condition
Psp_a_max=Ps_a_max*C(24)/ C(9);               % maximum output power happens at extreme
% conditions, the prop efficiency is considered to be the allowed
% worst, 60%. Also considers the average fuselage downwash for multicopters
n_a=(v_a_max+v_i_a)/(C(8)*C(10)*DsPm(4));   % rotation speed of the rotor, unit: rev/s
% J=v/(n*D), J advance ratio, v incoming flow velocity, D propeller diameter
% when operated at extreme conditions
Ps_max(3)=Psp_a_max;
n_max(3)=n_a;      % record candidate three for maximum motor output
% power and rotation speed
%% extreme condition analysis
% this condition is to work out the power and current requirement uncer the 
% extreme performance, it will then estimate the mass and price of motors, 
% ESCs and arms
[Ps_in_max, I]=max(Ps_max);   % find the max power and it's condition
n_max_req=n_max(I);   % rotation speed of the rotor, unit: rev/s

%Mach0.7 tip velocity check
omega_max=n_max_req*2*pi;
v_tip=omega_max*DsPm(4)/2;
if v_tip/C(4)>C(11)
    err_flag(9)=1;  % the prop tip velocity Mach number is exceeded, possible
    % noise and efficiency drop, should avoid. 
    
    err_percent=3;
    return;
end

kv_min=ceil(n_max_req*60*C(29)/(U*10))*10;       % kv is rpm per volt, so the rotation speed needs
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
[m_esc,pr_esc]=level3_esc(param_esc);
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
%% calculate the possible hovering and cruising endurance
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

Load(1)=Ts_max;              % the maximum thrust of a single rotor, unit N
Load(2)=max(n_max);       % the mini mum required rotation speed at maximum thrust, unit rev/s
Load(3)=max(Ps_max)/omega_max;     % the maximum torque, unit N.m=W/(rad/s)
Load(4)=m_mot;              % the motormass unit kg
Load(5)=n_he;               % hovering motor rotation speed
Load(6)=n_ce;               % cruising motor rotation speed
if Cstr(18)==1
    [m_arm, pr_arm,opt_geo]=level3_arm_propt(L_a);
else
    [m_arm,pr_arm, opt_geo]=level3_arm_mopt(L_a);
end
if DsPm(5)==2 && opt_geo(1)==1
    m_arm=m_arm*0.458;
end
%% final add up
M_tot=DsPm(1)*(m_mot+m_esc+m_p)+DsPm(6)*m_arm+m_b+m_body+m_av+m_LG;
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
if v_f_max<Cstr(1) 
    err_max_vf=1;   % the error flag of max forward speed not meet the mission requirment
    err_flag(3)=err_max_vf;
end
if v_a_max<Cstr(2) 
    err_max_va=1;   % error flag of max ascent speed not meet the mission requirement
    err_flag(4)=err_max_va;
end
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
%% mission constraints check
if v_f_max<max(v(:,1))
    missionerr(2)=1;   % the error flag of max forward speed not meet the mission requirment
end
if v_a_max<max(v(:,2))
    missionerr(3)=1;   % error flag of max ascent speed not meet the mission requirement
end
%% record the performance and design parameters of this multicopter

output=[Cstr(18),DsPm(1:3),DsPm(4)*1000,DsPm(5:7),DL,opt_geo(1),...
    opt_geo(2:5).*1000,opt_geo(6),L_d*1000,GTOM0*1000,M_tot*1000,...
    err_percent, Pr_tot,E_tot,FT_tot/60,FT_he/60,FT_ce/60,v_f_max,...
    v_a_max,m_mot*1000,m_esc*1000,m_p*1000,m_arm*1000,m_b*1000,m_body*1000,...
    m_av*1000,pr_mot,pr_esc,pr_p,pr_arm,pr_b,pr_body,P_av,Ps_in_max,...
    Psp_mot,omega_max,kv_min,I_dc_max,Cb,C_rate,H_lg_max*1000,err_flag,missionerr];

if sum(missionerr(2:3))~=0
    M_err=GTOM0*100;
end

end
function [Ts_g,n_g,Ps_g]=level3_general(var)
global C MidPara DsPm missionerr
% This function reads the velocity on two axis, and current GTOM, then
% estimate the required motor power, rotation speed and thrust.
% If it's the hovering condition, it will return the output power of the
% motor. However, if it's not hovering, it will only return the propeller
% output power, since the propeller efficiency would be different in
% differnt scenarios.

% there are three input variables: 1. forward speed, 2. ascending speed,
% 3. GTOM,
v_f=var(1);
v_a=var(2);
GTOW=var(3)*C(3);
% hovering condition
T_h=GTOW;      % equilibrium
Ts_h=T_h/DsPm(1);               % thrust of a single rotor
v_h=sqrt(Ts_h/(2*C(2,1)*MidPara(1,1)))  ;    % induced velocity in hover

if v_f==0 && v_a==0
    Ps_h=v_h*Ts_h;                  % required hovering power of a single rotor
    Psp_h=Ps_h*C(24)/C(6);            % prop input power in hover, fixed by FM
    % to consider the profile power, and fuselage downwash factor to
    % consider the frame area under the propeller.
    q=0;
    n_h=sqrt(Ts_h/(C(5,1)*C(2,1)*DsPm(4)^4));   % rotation speed of the rotor, unit: rev/s
else
    if v_f==0
        % vertical descent
        if v_a<0 && v_a/v_h>-2
           % fprintf("Warning: descent speed not appropriate,the multicopter will be in the VTR state, suggest move down-forward instead");
            missionerr(1)=1;
            %VTR; vortex ring state, it is an unstable state, the mulitcopter
            %can't get out from this state and will keep dropping
        end
        gama=atan(v_a/v_f);
        beta=0;         % tilt angle
        alpha=gama;    % prop angle of attack
        v_inf=v_a;      % velocity of comming wind
        q=C(2,1)*v_inf^2/2;       % dynamic pressure
    elseif v_a>=0
        gama=atan(v_a/v_f);     % forward flght, alpha==beta
        v_inf=sqrt(v_f^2+v_a^2);      % velocity of comming wind
        q=C(2,1)*v_inf^2/2;       % dynamic pressure
        alpha_ub=pi/2;      % max tilt angle
        alpha_lb=0;         % min tilt angle
        var1=[gama,v_inf,q,GTOW,alpha_ub,alpha_lb];
        alpha=level3_alpha(var1);  
    elseif v_a<0
        gama=atan(v_a/v_f);     % forward flght, alpha==beta
        v_inf=sqrt(v_f^2+v_a^2);      % velocity of comming wind
        q=C(2,1)*v_inf^2/2;       % dynamic pressure
        alpha_ub=pi/2;      % max tilt angle
        alpha_lb=-pi/2;         % min tilt angle
        var1=[gama,v_inf,q,GTOW,alpha_ub,alpha_lb];
        alpha=level3_alpha(var1);
    end
    D=q*C(26)*MidPara(2)*abs(sin(alpha)+cos(alpha)/C(27));  % Drag
    if pi/2-alpha<0.3
        Ts_g=(D+GTOW*sin(gama))/sin(alpha)/DsPm(1);
    else
    Ts_g=GTOW*cos(gama)/cos(alpha)/DsPm(1);
    end
%     T_g=Ts_g*DsPm(1);
    var_vi=[alpha, v_inf, Ts_g];
    v_i_g=level3_vi(var_vi);        % induced velocity    
    %P_prop=(D+GTOW*sin(gama))*v_inf; % =T_g*v_inf*sin(alpha), total power for movement
    Ps_prop=Ts_g*(v_i_g+v_inf*sin(alpha));    
    %Ps_prop=P_prop/DsPm(1);         % single rotor power for movement
    n_f=(v_inf*sin(alpha)+v_i_g)/(C(8)*C(10)*DsPm(4));   % rotation speed of the rotor, unit: rev/s
    
    if Ps_prop<0
        Ps_prop=0       % when descending up to a certain speed, the rotor gets
        % into wind mill brake state, the wind will work on the propeller to
        % deaccelerate the aircraft
    end
    
end




if v_f==0 && v_a==0
    P_mot=Psp_h;            % the motor output power in hovering
    n_g=n_h;
    Ts_g=Ts_h;
else
    P_mot=Ps_prop;      % the propelelr output power in all other conditions
    n_g=n_f;
end

Ps_g=P_mot;

end
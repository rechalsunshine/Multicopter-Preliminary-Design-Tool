function Ps_a=level3_ascent(v_a)
global C MidPara DsPm;
    % there are four input variables: 1. ascending speed, 
    % Middle Parameters: 1. single disk area--As, 
    % 2. top area--A1, 3. MTOW, 4, the hovering induced velocity
                         
    rho=C(2);                         % air density    
    q_a=rho*v_a^2/2;                   % dynamic pressure in ascending
    D=q_a*C(26)*MidPara(1,1);             % drag in ascending, q*Cd*A
    T_a=D+MidPara(3,1);                   % thrust in ascending, drag+GTOW
    Ts_a=T_a/DsPm(1);                % thrust of a single rotor, total thrust/N_rotor
%     if input_a(1)/input_a(6)>=0
%     v_i=-input_a(1)/2+sqrt((input_a(1)/2)^2+Ts_a/(2*C(2,1)*input_a(2)));   % the induced velocity 
%     % in ascending, it is derivated from T=2m'v_i=2(¦Ñ*A*v_in)*v_i,
%     % v_in=v_a+v_i, is the velocity of flow through the disk plane
%     elseif input_a(1)/input_a(6)<=-2
%         % when descending, the speed has to be large enough to avoid the
%         % vortex ring state
%         v_i=-input_a(1)/2-sqrt((input_a(1)/2)^2-Ts_a/(2*C(2)*input_a(2))); 
%     end
    Ps_a=Ts_a*v_a;   % required power for movement
%     P0_a=Ts_a*v_i;              % profile power
end
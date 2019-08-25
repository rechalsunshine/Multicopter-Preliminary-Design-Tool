function [m_body,pr_body,A_c, d_body, th_body]=level3_body(L_d)
global C DsPm MP MPr

% This function is used to calculate the mass and price of the central,
% which is modeled as three-layered round plate. two layer are used to
% fix all the arms, while the extra layer provide space for all the
% avionics and the battery.
% DsPm(5) is the current material in use: 1. Carbon Fiber, 2. Glass Fiber,
% 3. ABS Plastic, 4. Aluminium Alloy, 5. Wood
    d_body=0.2714*L_d^0.6718;
    th_body=MP(8,DsPm(5))*L_d*C(35);     % this parameter depends on materials
    % stronger material tend to need thinner plate. 
    A_c=pi*(d_body/2)^2 ;    % area of the central plate
    Vol=A_c*th_body;  % the total volumn of all the layers
    m_body=MP(1,DsPm(5))*Vol;  % call density from the material peroperty table
    pr_body=MPr(1,DsPm(5))*m_body;  % call unit price for plate material from bill table

end

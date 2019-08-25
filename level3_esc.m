function [m_esc,pr_esc,m_esc_chip]=level3_esc(param_esc)
    global C 

    % This function is used to calculate the mass and price of an ESC. The 
    % ESC is modelled as the chip plus signal wire, plus two input wire and
    % three output wire, plus any other radiator
     % There are two parameters used to estimate the mass and price of the 
    % ESC: 1. maximum input current, 2. the arm length
    
    %% old
    %d_in=0.0606*param_esc(1)+0.2941;    % a correlated relation to work out
    % the input wire diameter from the ampacity (continious maximum current)
    % unit: mm
    
    %n_in=floor(-39*log(d_in/0.127)/log(92)+36);     % get the id of gauge
    % from the estimated diameter. 
    
    % n_out=n_in+2;   % the out put wire is assumed to be one gauge smaller than
    % the input wire, since it's a three-phase current, the RMS current in 
    % each wire is smaller than the input current. And this is also what
    % all the ESC manufacturers do. 
    
    %m_wire_out=level3_wiremass(n_out,C(39));
    %% 
    n_in=floor(-4.016*log(param_esc(1))+29.998);
    %In order to meet the current
    % requirement, the gauge should be floor to the smaller integer (larger diameter) 
    
    % Assume use copper core for all the wires, and ignore the difference
    % between solid core and strand core. The density of copper is 8960 kg/m3. 
    rho_cp=C(37);
    m_wire_in=level3_wiremass(n_in, param_esc(2));
    m_esc_chip=0.8893*param_esc(1)/1000;  % work out the chip mass from correlation
    m_esc=2*m_wire_in+m_esc_chip;  % two input wire, three output wire
    % plus the chip
    
    pr_esc=8.6812*exp(0.0166*param_esc(1));

end

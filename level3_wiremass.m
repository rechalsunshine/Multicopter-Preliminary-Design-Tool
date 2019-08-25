
    function m_wire=level3_wiremass(n_wire, l_wire)
    global C
    
    % This function takes two input variables: 1. the wire gauge ID, 2. the
    % length of the wire. The output is the mass of this piece of wire.
    d_wire=0.127*92^((36-n_wire)/39);
    r=1.7086*d_wire^(-0.346);   % to work out the (core+insulation)/core ratio
    % of the given diameter. This coefficient is used to estimate the full
    % wire mass from only the density of the core. 
    A_wire=pi*(d_wire/2000)^2;   % calculate the wire cross section area
    rho_l=C(37)*A_wire*r;        % calculate the length density. unit: kg/m
    m_wire=rho_l*l_wire;        % calculate the mass for a given length of wire
    end
function [m_b,pr_b]=level3_bat(E_b_req)
    

    % This function is used to calculate the mass and price of a battery.
         
    m_b=0.0059*E_b_req;
    pr_b=0.6337*E_b_req;

end

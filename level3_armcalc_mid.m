%%% This function is a wrap of the LaWorkout function, it only output the
%%% mass of the arm, and this function is used for the mass optimisation.

function m_a=level3_armcalc_mid(property)


Arm_output=level3_armcalc(property);

if Arm_output.frequency_error_x<0 && Arm_output.frequency_error_y<0 &&...
        Arm_output.deflection_error<1.2e-6 && Arm_output.shear_stress<0 &&...
        Arm_output.normal_compression_stress<0 && Arm_output.normal_tension_stress<0
     
        m_a=Arm_output.mass;
else
    m_a=Arm_output.mass*100;
end


end
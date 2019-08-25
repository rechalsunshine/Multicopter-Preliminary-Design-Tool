%%% This function is a wrap of the LaWorkout function, it only output the
%%% mass of the arm, and this function is used for the mass optimisation.

function m_a=level3_armcalc_fmin(property)

global input

for i=2:5
    input(i)=property(i-1);  %test
end


Arm_output=level3_armcalc(input);

if Arm_output.frequency_error_x<0 && Arm_output.frequency_error_y<0 &&...
        Arm_output.deflection_error<1.2e-6 && Arm_output.shear_stress<0 &&...
        Arm_output.normal_compression_stress<0 && Arm_output.normal_tension_stress<0
     
        m_a=Arm_output.mass;
else
    m_a=Arm_output.mass*1000;
end


end
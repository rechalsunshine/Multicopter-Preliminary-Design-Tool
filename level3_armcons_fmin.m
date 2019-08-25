%%% This function is a wrap of the LaWorkout function, it outputs the
%%% constrainsts.

function [cons,ceq]=level3_armcons_fmin(property)
global input
for i=2:5
    input(i)=property(i-1);
end

Arm_output=level3_armcalc(input);
cons(1)=Arm_output.frequency_error_x;      % Consider the first order frequency on both x and y axis
cons(2)=Arm_output.frequency_error_y;      % the motor frequency should be smaller than these frequencies 
cons(3)=Arm_output.deflection_error;        % deflection constraint
cons(4)=Arm_output.shear_stress;
cons(5)=Arm_output.normal_compression_stress;
cons(6)=Arm_output.normal_tension_stress;

ceq=[];

end
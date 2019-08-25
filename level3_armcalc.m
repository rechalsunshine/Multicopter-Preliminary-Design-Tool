function Arm=level3_armcalc(property)

global C MP Load err_flag
Arm = struct('mass', 1,'price', 1,'frequency_error_x', 1,'frequency_error_y', 1,'deflection_error', 1,...
    'shear_stress', 1,'normal_compression_stress', 1,'normal_tension_stress', 1);


%read fixed parameters from a given scenario
%FixedParameters=xlsread('FixedParameters.xlsx',MissionName);
% read fixed parameters from a given scenario

%%%% read geometry parameters

SFL=property(1);         % Shape flag 1:rectangular plate, 2: hollow rectangular, 3: hollow oval, 4: T section, 5: Channel section, 6: I section
w=property(2);              % width,unit m
h=property(3);             % height,unit m
t=property(4);          % thickness,unit m
L=property(5);          % Length of the amr, unit m
MFL=property(6);        % Material Flag 1: Carbon Fibre, 2: Glass Fibre, 3: ABS plastic, 4: Aluminium

%read material parameters

rho=MP(1,MFL);            % material density, kg/m3
E=MP(2,MFL);              % Young's modulus, Pa
G=MP(3,MFL);              % shear modulus, Pa
xigma_tn=MP(4,MFL);       % Tension yield/ultimate strength, Pa
xigma_cp= MP(5,MFL);      % Compression yield/ultimate strength, Pa
tao=MP(6,MFL);            % Shear yield/ultimate strength, Pa
Pr_uni=MP(7,MFL);         %??? unit price of the material

if t>w/2 || t>h/2
    m_arm=5;
    Pr_arm=1000;
    fex=1;
    fey=1;
    de=1;
    tao_err=1e8;
    xigma_err=[1e8,1e8];
else
        %% %%%% calculation of the load %%%%%%%%
    T_max=Load(1);         % unit N
    n_max=Load(2);         % unit rev/s=Hz, workout the rotation speed
    Q_max=Load(3);         %unit N.m, workout the torque of the motor
    
    Fx=T_max-Load(4)*C(3);
    Mx=Fx*L;
    My=Q_max;    

    %% %%%%calculation of geometry properties, including cross section area (A) and area moment of inertia (Ix) %%%%
    % syms w h t   for test
    switch SFL
        case 1  % rectangular          
            Ix=w*h^3/12;        % x axis area moment of inertia
            Iy=h*w^3/12;        % y axis area moment of inertia
            Ixy=0;
            A=w*h;              % cross section area
            Sx=w*h/2*h/4;       % the area moment to the dangerous area (for shear stress calculation)
            we=w;
            Ae=Ix*we/Sx;         % equivalent area for shear stress, The explanation is below
            h_mc=h/2;
            h_mt=-h/2;
            w_mc=w/2;
            w_mt=-w/2;
        case 2  % hollow rectangular (box)
            Ix=(w*h^3-(w-2*t)*(h-2*t)^3)/12;
            Iy=(h*w^3-(h-2*t)*(w-2*t)^3)/12;
            Ixy=0;
            A=w*h-(w-2*t)*(h-2*t);
            Ae=2/3*(w*h^3-(w-2*t)*(h-2*t)^3)/(h^2+2*(w-2*t)*(h-t));
            h_mc=h/2;
            h_mt=-h/2;
            w_mc=w/2;
            w_mt=-w/2;
        case 3  % hollow oval (oval cilindar)
            Ix=pi*(w*h^3-(w-2*t)*(h-2*t)^3)/64;
            Iy=pi*(h*w^3-(h-2*t)*(w-2*t)^3)/64;
            Ixy=0;
            A=pi*(w*h-(w-2*t)*(h-2*t))/4;
            Ae=3*pi/4*((w/2)*(h/2)^3-(w/2-t)*(h/2-t)^3)*t/((w/2)*(h/2)^2-(w/2-t)*(h/2-t)^2);  
            k=-(My*Ix+Mx*Ixy)/(Mx*Iy+My*Ixy);
%              syms w_m h_m;       
%             eq1=k==-h^2*w_m/(w^2*h_m);
%             eq2=w_m^2/(w/2)^2+h_m^2/(h/2)^2==1;
% %           eq1=k+h^2*w_mc/(w^2*h_mc);
% %           eq2=w_mc^2/(w/2)^2+h_mc^2/(h/2)^2==1;
%             [h_m,w_m]=solve(eq1,eq2);       % solve a two parameter power two 
%             % equations to find the coordinate of the tangent point on the oval
%             % There are two point where max tension stress and compression
%             % stress happens
            w_mt=w^2*k/(2*sqrt(w^2*k^2+h^2));
            h_mt=-h^2/(2*sqrt(w^2*k^2+h^2));
            w_mc=-w^2*k/(2*sqrt(w^2*k^2+h^2));
            h_mc=h^2/(2*sqrt(w^2*k^2+h^2));
        case 4  % "T" section
            e1=(w*t+h^2-t^2)/(h-t+w)/2;
            Ix=(w*e1^3-(w-t)*(e1-t)^3+t*(h-e1)^3)/3;
            Iy=(t*w^3+(h-t)*t^3)/12;
            Ixy=0;
            A=w*t+(h-t)*t;      
            Ae=2/3*(w*e1^3-(w-t)*(e1-t)^3+t*(h-e1)^3)/(h-e1)^2;
            w_mt=-t/2;
            h_mt=e1-h;
            w_mc=t/2;
            h_mc=e1;
        case 5  % Channel section
            e1=(w*t+2*(h^2-t^2))/(2*h+w-2*t)/2;
            Ix=(w*e1^3-(w-2*t)*(e1-t)^3+2*t*(h-e1)^3)/3;
            Iy=(h*w^3-(h-t)*(w-2*t)^3)/12;
            Ixy=0;
            A=w*t+(h-t)*2*t;
            Ae=2/3*(w*e1^3-(w-2*t)*(e1-t)^3+2*t*(h-e1)^3)/(h-e1)^2;
            w_mt=-w/2;
            h_mt=e1-h;
            w_mc=w/2;
            h_mc=e1;
        case 6 % "L" section
            ey=(w*t+h^2-t^2)/(h-t+w)/2;
            ex=(h*t+w^2-t^2)/(h-t+w)/2;
            Ix=(w*ey^3-(w-t)*(ey-t)^3+t*(h-ey)^3)/3;            
            Iy=(h*ex^3-(h-t)*(ex-t)^3+t*(w-ex)^3)/3;
            Ixy=w*h*t*(w-t)*(h-t)/4;
            A=w*t+(h-t)*t;      
            Ae=2/3*(w*ey^3-(w-t)*(ey-t)^3+t*(h-ey)^3)/(h-ey)^2;
             
            w_mt=ex-t;      % assume torque is to the right ->, and thrust is to the up,then 
            % the combined torque is to the upright. the optimal shape of
            % the beam cross section is like "7" or "L", not the reverse.
            h_mt=ey-h;
            w_mc=ex;
            h_mc=ey;
%             w_mt=-ex;     % This is when place the crosssection like "J"
%             or "F".assume the same load direction
%             h_mt=ey-h;
%             w_mc=w-ex;
%             h_mc=ey;
       
        case 7   % "I" section
            Ix=(w*h^3-(w-t)*(h-2*t)^3)/12;
            Iy=(2*t*w^3+(h-2*t)*t^3)/12;
            Ixy=0;
            A=w*h-(w-t)*(h-2*t);
            Ae=2/3*(w*h^3-(w-t)*(h-2*t)^3)/((h-2*t)^2+4*w*(h-t));
            h_mc=h/2;
            h_mt=-h/2;
            w_mc=w/2;
            w_mt=-w/2;

        otherwise % "H" section
            Iy=(h*w^3-(h-t)*(w-2*t)^3)/12;
            Ix=(2*t*h^3+(w-2*t)*t^3)/12;
            Ixy=0;
            A=w*h-(h-t)*(w-2*t);
            Ae(1)=2/3*w*(2*h^3+w*t^2-2*t^3)/(2*h^2+w*t-2*t^2);      % middle plane
            Ae(2)=2/3*t*(2*h^3+w*t^2-2*t^3)/(h-t)^2;                % edge of middle plane
            h_mc=h/2;
            h_mt=-h/2;
            w_mc=w/2;
            w_mt=-w/2;
    end
    
    %% %%%% calculation of the beam performance (constraints) %%%%%%%
    fnx=sqrt((3*E*Ix/L^3)/(rho*A*L*33/140+Load(4)))/(2*pi);        %???? the equation should be corrected for a cantilever beam with tipped mass at the free end
    fny =sqrt((3*E*Iy/L^3)/(rho*A*L*33/140+Load(4)))/(2*pi);
    FEX(1)=C(32)-abs(fnx-n_max)/fnx; % ????the frequency error from the natural frequency should be larger than 0.3
    FEX(2)=C(32)-(fnx-Load(5))/fnx;
    FEX(3)=C(32)-(fnx-Load(6))/fnx;
    FEY(1)=C(32)-abs(fny-n_max)/fny;
    FEY(2)=C(32)-(fny-Load(5))/fny;
    FEY(3)=C(32)-(fny-Load(6))/fny;
    
    fex=max(FEX);
    fey=max(FEY);
    if fnx<n_max || fny<n_max
        err_flag(8)=1;  % the natural frequency of the beam is smaller than
        % the maximum rotation speed, 
    end
    
    dmax_x=Fx*L^3/(3*E*Ix);       % the maximum x deflection of the beam
    dmax_y=My*L^2/(2*E*Iy);        % the maximum y deflection of the beam
    dmax=sqrt(dmax_x^2+dmax_y^2);
    d_limit=L/C(33);                  % ???? work out the deflection limit
    de=dmax-d_limit;              % the maximum deflection error
        
    for i=1:length(Ae)
        tao_dan(i)=Fx/Ae(i);    % calculate shear stress on all the dangerous planes
    end
    tao_max=max(tao_dan);                   % the maximum flexural shear stress caused by thrust
    tao_err=tao_max-tao/C(34);         % error between allowed maximum shear stress vs the max shear stress on the cross section
    % tao=Fx*Sx/(Ix*we), [Fx] is the shear force, [Sx] is the area moment to
    % the dangerous plane, [Ix] is the area moment of inertia, and [we] is the
    % equivalent width. This canbe written as tao=Fx/Ae, so for each shape,
    % the equivalent stress area Ae=Ix*we/Sx.
    
    xigma_max_com=abs(((Mx*Iy+My*Ixy)*h_mc-(My*Ix+Mx*Ixy)*w_mc)/(Ix*Iy-Ixy^2));
    xigma_max_ten=abs(((Mx*Iy+My*Ixy)*h_mt-(My*Ix+Mx*Ixy)*w_mt)/(Ix*Iy-Ixy^2));
    %xigma_max=Mx*h_m/Ix-My*w_m/Iy;    % the maximum stress in the beam, it is 
    % the overall stress under both thrust and torque applied on the beam.
    % In this euqation, [h_m] and [w_m] are the distances from the
    % dangerous point (max stress point) to the centroid point. 
    
    % xigma_maxq=Q_max*(w/2)/Iy;   %??? the maximum stress caused by torque
    xigma_err(1)=xigma_max_ten-xigma_tn/C(34); % ????the error between the allowed normal stress vs the actual maximum normal stress
    xigma_err(2)=xigma_max_com-xigma_cp/C(34);
    
    
    m_arm=rho*A*L;          % final mass of the arm
    Pr_arm=Pr_uni*m_arm;      % final price of the arm
end

%% %%%%final output of the function %%%%

Arm = struct('mass', m_arm,'price', Pr_arm,'frequency_error_x', fex,'frequency_error_y', fey,'deflection_error', de,...
    'shear_stress', tao_err,'normal_compression_stress', xigma_err(2),'normal_tension_stress', xigma_err(1));
end
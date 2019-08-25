function a=level3_alpha(var)
global C MidPara

% four inputs:1. gama, 2. v_inf, 3. q, 4. GTOW, 5. alpha_ub, 6, alpha_lb
if var(1)==pi/2 || var(1)==-pi/2
    alpha=var(1);
else
Cons=var(3)*C(26)*MidPara(2)/(var(4)*cos(var(1)));

alpha =var(5)*0.999 ;   % the start point of alpha
alpha_old = 100;    % the alpha value in the last iteration
iter = 0;   % initialise the iteration flag
while abs((alpha_old-alpha)/alpha) > 10^-9 && alpha ~= 0 || alpha>var(5) || alpha < var(6)
    alpha_old = alpha;
    f_A=(sin(alpha)+cos(alpha)/C(27));
    f_Ad=(cos(alpha)-sin(alpha)/C(27));
    if f_A>0
    alpha=alpha-(Cons*f_A-tan(alpha)+tan(var(1)))/(Cons*f_Ad-1/cos(alpha)^2);
    else
        alpha=alpha-(-Cons*f_A-tan(alpha)+tan(var(1)))/(-Cons*f_Ad-1/cos(alpha)^2);
    end
    % above is the Newton-Raphson method, alpha_n+1=alpha_n-(f(alpha)/f'(alpha))_n;
    % f(alpha)=Cons*(sin(alpha)+cos(alpha)/n)+tan(gama)-tan(alpha)==0. 
    % f'(alpha)=Cons*(cos(alpha)-sin(alpha)/n)-1/cos(alpha)^2
    % C=Cd*A1/MTOW
    iter = iter + 1;
    alpha_d=alpha*180/pi;
%     fprintf('Iteration %d: alpha=%.20f degree, %.20f rad, err=%.20f\n', iter, alpha_d,alpha, alpha_old-alpha);
    %pause;
    if alpha>pi/2 
        alpha=pi-alpha;
    elseif alpha<-pi/2
        alpha=pi+alpha;
    end
end
gama_d=var(1)*180/pi;
%     fprintf('gamma= %.20f, beta=%.20f', gama_d,alpha_d-gama_d );
end
a=alpha;

end


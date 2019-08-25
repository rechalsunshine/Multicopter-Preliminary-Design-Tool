function v_i=level3_vi(var_vi)
global C MidPara

% four inputs:1. alpha, 2. v_inf, 3. Ts, 

vi0 =20 ;   % the start point of alpha
vi=vi0;
vi_old = 100;    % the alpha value in the last iteration
iter = 0;   % initialise the iteration flag
while abs((vi_old-vi)/vi) > 10^-9 && vi ~= 0 
    vi_old = vi;
    vi=vi-(vi^4+2*var_vi(2)*vi^3*sin(var_vi(1))+var_vi(2)^2*vi^2-(var_vi(3)/(2*C(2)*MidPara(1)))^2)/...
        (4*vi^3+6*var_vi(2)*sin(var_vi(1))*vi^2+2*var_vi(2)^2*vi);
%     fprintf('Iteration %d: vi=%.20f m/s, err=%.20f\n', iter, vi, vi_old-vi);
    % above is the Newton-Raphson method, alpha_n+1=alpha_n-(f(alpha)/f'(alpha))_n;
    % f(x)=x^4+2*v_inf(i)*x^3*sin(alpha(i))+v_inf(i)^2*x^2-(T(i)/C)^2==0. 
    % f'(x)=4x^2+6v_inf*sin(alpha)*x^2+2v_inf^2*x
    % C=Cd*A1/MTOW
    iter = iter + 1;
%     if vi<0
%         vi=vi0*iter*2;
%     end  
end
v_i=vi;

end


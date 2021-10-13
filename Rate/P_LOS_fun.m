function [P_LOS] = P_LOS_fun(t,omega,phi,lambda_b,E_L,L_min,L_max)

a = t.*abs(sin(phi-omega))./sin(phi);

[M, N] = size(omega);

for m = 1:M
    for n = 1:N
        if 0 <= a(m,n) && a(m,n) < L_min
            P_LOS(m,n) = t(m,n)*exp(-lambda_b*t(m,n)^2*abs(sin(phi-omega(m,n)))*sin(omega(m,n))/(2*sin(phi)));
        elseif (L_min <= a(m,n) && a(m,n) < L_max)
            P_LOS(m,n) = t(m,n)*exp(-lambda_b/(L_max-L_min)*((a(m,n)^2-L_min^2)/2*t(m,n)*sin(omega(m,n))-(a(m,n)^3-L_min^3)/6*sin(phi)*sin(omega(m,n))/abs(sin(phi-omega(m,n)))+(L_max-a(m,n))*t(m,n)^2*abs(sin(phi-omega(m,n)))*sin(omega(m,n))/(2*sin(phi))));
        elseif (L_max <= a(m,n))
            P_LOS(m,n) = t(m,n)*exp(-lambda_b*(E_L*t(m,n)*sin(omega(m,n))-1/(L_max-L_min)*(L_max^3-L_min^3)/6*sin(phi)*sin(omega(m,n))/abs(sin(phi-omega(m,n)))));
        end 
%         P_LOS(m,n) = t(m,n)*exp(-lambda_b*t(m,n)^2*abs(sin(phi-omega(m,n)))*sin(omega(m,n))/(2*sin(phi)))*(0 <= a(m,n) && a(m,n) < L_min)...
%                      +t(m,n)*exp(-lambda_b/(L_max-L_min)*((a(m,n)^2-L_min^2)/2*t(m,n)*sin(omega(m,n))-(a(m,n)^3-L_min^3)/6*sin(phi)*sin(omega(m,n))/abs(sin(phi-omega(m,n)))+(L_max-a(m,n))*t(m,n)^2*abs(sin(phi-omega(m,n)))*sin(omega(m,n))/(2*sin(phi))))*(L_min <= a(m,n) && a(m,n) < L_max)...
%                      +t(m,n)*exp(-lambda_b*(E_L*t(m,n)*sin(omega(m,n))-1/(L_max-L_min)*(L_max^3-L_min^3)/6*sin(phi)*sin(omega(m,n))/abs(sin(phi-omega(m,n)))))*(L_max <= a(m,n));
%         
    end
end
end


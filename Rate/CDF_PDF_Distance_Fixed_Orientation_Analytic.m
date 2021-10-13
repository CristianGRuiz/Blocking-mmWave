clc
close all;
clear all;

% ACCESS POINTS
% The density of the access points
lambda_AP = 1e-4;

% BLOCKING ELEMENTS
% The line segments model is analyzed. The elements shapes are
% characterized by their lengths L with L~U[L_min,L_max]
[E_L, lambda_b, NumberBuildings, A_T, F_L_emp, l_emp_cdf] = obtain_average_densities_circumscript();

L_min = 0.1;
L_max = 57.01;

% POWER PARAMETERS.
P_U = 33;                           % Tx Power of the User in dBm.
    P_U = 10^(P_U/10)/1000;         % Tx Power of the User in Watts.
BW = 10e6;                          % Bandwidth of the system.
S_0 = -174;                         % Noise Power Spectral Density in 
                                    % dBm/Hz.
    S_0 = 10^(S_0/10)/1000;         % Noise Power Spectral Density in 
                                    % Watts/Hz.
sigma_2 = S_0*BW;                   % Noise Power in Watts.
L_ref = 25.6;                       % Path Loss at Reference distance in 
                                    % dB.
    L_ref = 10^(L_ref/10);          % Path Loss at Reference distance in 
                                    % linear scale.
alpha = 4;                          % Path Loss Exponent.
gamma = P_U/(sigma_2*L_ref);        % SNR at the reference distance.
chi = 0.577;                        % Euler-Mascheroni constant.
rho = exp(-chi);

% We define some parameters for the bounds and approximations for the sinus
m = (96*pi-24)/(4*pi^4-3*pi^2);
n = (8-m*pi^2)/(4*pi);

% We create an array d for evaluating the expression
d = 0:10:1000;

% We obtain the expressions for F_D_ub and F_D_lb
warning('off','all')
% We obtain the CDF for the scenario without blockages (NB)
F_D_NB = 1 - exp(-lambda_AP*pi*d.^2);
for p=1:length(d)
    F_D_ub(p) = 1 - exp(-4*lambda_AP*integral2(@(x,phi) x.*exp(-lambda_b*E_L*x.*sin(phi)),0,d(p), 0, pi/2, 'RelTol', 1e-1, 'AbsTol', 1e-1));
    F_D_app(p) = 4*lambda_AP*integral2(@(x,phi) g_x_phi_fun(x,phi,lambda_AP,lambda_b,E_L,L_min,L_max), 0, d(p), 0, pi/2, 'RelTol', 1e-1, 'AbsTol', 1e-1);
end
warning('on','all')
F_D_ub_app = 1 - exp(-4*lambda_AP/((lambda_b*E_L)^2*m*n*(n+m*pi/2))*(m*pi/2+n*exp(-lambda_b*E_L*d*(n+m*pi/2))-(n+m*pi/2)*exp(-lambda_b*E_L*d*n)));

% Now, we obtain the results for the statistics of the lower bound of the
% rate.
r = 0:25;
warning('off','all')
for p=1:length(r)
    F_R_lb_app(p) = 1-4*lambda_AP*integral2(@(x,phi) g_x_phi_fun(x,phi,lambda_AP,lambda_b,E_L,L_min,L_max), 0, ((gamma*rho)/(exp(r(p))-1))^(1/alpha), 0, pi/2, 'RelTol', 1e-3, 'AbsTol', 1e-3);
end

save('CDFs_PDFs_Distance_Fixed_Orientation_Analytic.mat','d','r','F_D_NB','F_D_ub','F_D_ub_app','F_D_app','F_R_lb_app');
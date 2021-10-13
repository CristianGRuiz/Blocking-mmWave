clc
close all;
clear all;

% In this script, for the case where the orientation of the blocking
% elements is fixed, we obtain the distance to the closest visible AP to
% the reference user which is placed at the origin for simplicity purposes.
% Then, the empirical CDF and PDF of this random variable are derived.

% SIMULATION'S PARAMETERS
N_run = 10000;            % Total number of simulations.
D_max = 3000;            % Maximum considered distance.

% AP DEPLOYMENT'S PARAMETERS
lambda_AP = 1e-4;       % Density of APs.

% BLOCKINGS' PARAMETERES
% The line segments model is analyzed. The elements shapes are
% characterized by their lengths L with L~U[L_min,L_max].
[E_L, lambda_b, NumberBuildings, A_T, F_L_emp, l_emp_cdf] = obtain_average_densities_circumscript();

L_min = 0.1;
L_max = 57.01;

% Initialization of all the variables.
D_emp = Inf*ones(1,N_run);
R_emp = zeros(1,N_run);

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

% Iterate through each simulation.
for n_run = 1:N_run
%     close all
%     clc
    % ACCESS POINTS
    % Calculate both the number and the characteristics of the APs.
    % Calculate the number of APs of every simulation, according to fact
    % that N_AP is a Poisson variable, and its mean value can be obtained
    % by multplying the density for the area.
    N_AP = poissrnd(lambda_AP*pi*D_max^2);
    % Distance from the different APs to the origin, in ascending order.
    d_AP = sort(D_max*sqrt(rand(1,N_AP)));
    % Angle of the different APs with respect to the horizontal.
    phi_AP = 2*pi*rand(1,N_AP);
    % Coordinates of the center of the APs.
    x_AP = d_AP.*cos(phi_AP);        y_AP = d_AP.*sin(phi_AP);
    
    % BLOCKING ELEMENTS
    % Calculate both the number and the characteristics of the blocking
    % elements.
    % Calculate the number of buildings of every simulation, according to
    % the fact that B is a Poisson variable, and its mean value can be
    % obtained by multplying the density for the area. A margin of L_max/2
    % is given to ensure that the density of blockings is the same in the
    % borders.
    B = poissrnd(lambda_b*pi*(D_max + L_max/2)^2);
    % Distance from the center of the different blockings to the origin,
    % again with the same margin of L_max/2, in ascending order.
    d_b = sort((D_max + L_max/2)*sqrt(rand(1,B)));
    % Angle of the different blockings with respect to the horizontal.
    phi_b = 2*pi*rand(1,B);
    % Length of every building.
    l = L_min + (L_max - L_min)*rand(1,B);
    % Coordinates of the center of the blocking elements.
    x_b = d_b.*cos(phi_b);        y_b = d_b.*sin(phi_b);
    
%     figure(1)
%     hold on;
%     plot(0,0,'ro','MarkerSize',10,'LineWidth',2);
%     text(0,0,'Reference User')
%     for b = 1:B
%         plot([x_b(b) - l(b)/2, x_b(b) + l(b)/2],[y_b(b) y_b(b)],'k-','LineWidth',2);
%     end
%     for n_AP = 1:N_AP
%         plot(x_AP(n_AP),y_AP(n_AP),'bo','MarkerSize',10,'LineWidth',2);
%     end
%     xlim([-D_max D_max]);      ylim([-D_max D_max]);
%     daspect([1 1 1])
    
    % Iterate through each AP until it is visible to the reference user.
    for n_AP = 1:N_AP
        % Obtain the intersection point between the equations that
        % represent the link from the reference user to the n_AP-th AP and
        % all the blocking elements.
        x = x_AP(n_AP)*y_b/y_AP(n_AP);      y = y_b;
        
        % This point must meet 3 requirements:
        % 1: This point is contained by the building (that is, the point is
        % closer than l/2 to the center of the blocking element.
        aa1 = abs((x+1i*y)-(x_b+1i*y_b)) <= l/2;
        % 2: The following 2 conditions ensure that this point is between
        % the origin and the AP.
        aa2 = abs((x+1i*y)) <= abs(x_AP(n_AP)+1i*y_AP(n_AP));
        aa3 = abs((x+1i*y)-(x_AP(n_AP)+1i*y_AP(n_AP))) <= abs(x_AP(n_AP)+1i*y_AP(n_AP));
        
        % If there is at least one blocking element that blocks the link,
        % the intersection point (x,y) meets the three conditions. Then, if
        % and only if these conditions at the same time are not satisfied
        % by any point, the AP is visible to the origin.
        if ~sum(aa1.*aa2.*aa3)
            D_emp(n_run) = d_AP(n_AP);
            h = 1/sqrt(2)*(randn(1) + 1i*randn(1));
            R_emp(n_run) = log(1+(D_emp(n_run))^(-alpha)*abs(h)^2*gamma);
%             d_AP(n_AP)
%             plot(x_AP(n_AP),y_AP(n_AP),'go','MarkerSize',10,'LineWidth',2);
%             plot([0,x_AP(n_AP)],[0,y_AP(n_AP)],'g--','LineWidth',2);
            break
        end
    end  
end

% Obtain the CDF
[F_D_emp,d_emp_cdf] = ecdf(D_emp);
[F_R_emp,r_emp_cdf] = ecdf(R_emp);

% We rotate the resulting vectors.
F_D_emp = F_D_emp.';
d_emp_cdf = d_emp_cdf.';
F_R_emp = F_R_emp.';
r_emp_cdf = r_emp_cdf.';

% Obtain the PDF
[f_D_emp,d_emp_pdf] = ksdensity(D_emp(D_emp ~= Inf));
[f_R_emp,r_emp_pdf] = ksdensity(R_emp);

% Save the results
save('CDFs_PDFs_Distance_Fixed_Orientation_Simulation.mat','F_D_emp','d_emp_cdf','F_R_emp','r_emp_cdf','f_D_emp','d_emp_pdf','f_R_emp','r_emp_pdf');
clc
close all;
clear all;

load('CDFs_PDFs_Distance_Fixed_Orientation_Analytic.mat');
load('CDFs_PDFs_Distance_Fixed_Orientation_Simulation.mat');

% Plot section
lw = 1;             % Line Width
fs = 15;            % Font Size
ms = 12;            % Marker Size
md = 50;            % Marker Distance

figure(1)
set(gca,'defaultTextInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')
hold on
plot(d,F_D_NB,'b-','LineWidth',lw,'HandleVisibility','off');  plot(d(1:md:length(d)),F_D_NB(1:md:length(d)),'bs','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'bs-','LineWidth',lw,'MarkerSize',ms);
plot(d,F_D_ub,'g-','LineWidth',lw,'HandleVisibility','off');  plot(d(round(2*md/5):md:length(d)),F_D_ub(round(2*md/5):md:length(d)),'g>','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'g>-','LineWidth',lw,'MarkerSize',ms);
plot(d,F_D_ub_app,'k-','LineWidth',lw,'HandleVisibility','off');  plot(d(round(3*md/5):md:length(d)),F_D_ub_app(round(3*md/5):md:length(d)),'ko','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'ko-','LineWidth',lw,'MarkerSize',ms);
plot(d,F_D_app,'m-','LineWidth',lw,'HandleVisibility','off');  plot(d(round(4*md/5):md:length(d)),F_D_app(round(4*md/5):md:length(d)),'m+','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'m+-','LineWidth',lw,'MarkerSize',ms);
plot(d_emp_cdf,F_D_emp,'r--','LineWidth',2);
legend('$F_D(d)$ without blockages',...
       '$\overline{F_D(d)}$ with $P(Z | T_1)$',... 
       'Appr. of $\overline{F_D(d)}$ with appr. on $P(Z | T_1)$',...       
       'Appr. of $F_D(d)$ with $g(x,\phi)$',...
       'Empiric $F_D(d)$',...
       'Interpreter','latex','location','southeast')
title('Bounds and approximations for $F_D(d)$')
xlabel('Distance from the RU at O to the closest visible AP, $d$[m]')
ylabel('$F_D(d)$')
set(gca,'fontsize',fs)
xlim([0 min([max(d), max(d_emp_cdf)])]);
ylim([0 1]);
grid on
grid minor
hold off;

md = 5;            % Marker Distance

figure(2)
set(gca,'defaultTextInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')
hold on
plot(d,F_D_NB,'b-','LineWidth',lw,'HandleVisibility','off');  plot(d(1:md:length(d)),F_D_NB(1:md:length(d)),'bs','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'bs-','LineWidth',lw,'MarkerSize',ms);
plot(d,F_D_ub,'g-','LineWidth',lw,'HandleVisibility','off');  plot(d(round(2*md/5):md:length(d)),F_D_ub(round(2*md/5):md:length(d)),'g>','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'g>-','LineWidth',lw,'MarkerSize',ms);
plot(d,F_D_ub_app,'k-','LineWidth',lw,'HandleVisibility','off');  plot(d(round(3*md/5):md:length(d)),F_D_ub_app(round(3*md/5):md:length(d)),'ko','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'ko-','LineWidth',lw,'MarkerSize',ms);
plot(d,F_D_app,'m-','LineWidth',lw,'HandleVisibility','off');  plot(d(round(4*md/5):md:length(d)),F_D_app(round(4*md/5):md:length(d)),'m+','LineWidth',lw,'MarkerSize',ms,'HandleVisibility','off'); plot(nan,nan,'m+-','LineWidth',lw,'MarkerSize',ms);
plot(d_emp_cdf,F_D_emp,'r--','LineWidth',2);
set(gca,'fontsize',fs)
xlim([0 150]);
ylim([0 0.7]);
grid on
grid minor
hold off; 

md = 50;            % Marker Distance

figure(3)
set(gca,'defaultTextInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')
hold on
plot(r,F_R_lb_app,'b-','LineWidth',2);
plot(r_emp_cdf,F_R_emp,'r--','LineWidth',2);
legend('Appr. of $F_{\underline{R}}(r)$ with $g(x,\phi)$',...
       'Empiric $F_R(r)$',...
       'Interpreter','latex','location','southeast')
title('CDF of the rate')
xlabel('Rate [nats/Hz]')
ylabel('$F_R(r)$')   
set(gca,'fontsize',fs)
xlim([0 min([max(r), max(r_emp_cdf)])]);
ylim([0 1]);
grid on
grid minor
hold off;      
%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 2, Problem 1d
% September 17, 2026

%--------------------------------------------------------------------------
%% Plotting speed against threshold in a model for a cable with idealized excitable nonlinearity
D=1;


c_func = @(a, D) sqrt(D)*(1-2.*a)./(sqrt(a.*(1-a)));

a_span = 0.0001:0.000001:1-0.0001;
c_vals = c_func(a_span, D);
[~, index] = min(abs(c_vals - 2));
cross_a = a_span(index);


figure (1)
plot(a_span,c_vals, 'LineWidth', 1)
x = xlim;
y = ylim;
hold on
%plot(t,u2, 'LineWidth', 1)
%plot(t,u3, 'LineWidth', 1)
%y = ylim; % current y-axis limits
%plot([tpeak1 tpeak1],[y(1) y(2)], 'LineWidth',1) % visual check
plot([x(1) x(2)], [2 2], 'k--', 'LineWidth',1)
plot([cross_a cross_a], [y(1) y(2)], 'r--', 'LineWidth',1)
hold off
title('Propogation Speed vs. Excitation Threshold, $D=1$, min(a)=0.0001', 'interpreter','latex','FontSize', 14)
legend('c', 'passive apparent speed', 'crossing value')
xlabel('$a$ (excitation threshold)','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$c$ (propogation speed)','interpreter','latex', 'FontWeight','bold','FontSize',16)




%==========================================================================








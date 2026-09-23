%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 2, Problem 2c
% September 17, 2026

%--------------------------------------------------------------------------
%% Finding the voltage V for which the slope dIN/dV vanishes
Mg = 1;
EN = 0;
V0 = 16.13;
g=1;


dIN_func = @(V) g*(1+(Mg/3.57).*exp(-V/V0)).^(-1).*(1+((V-EN).*(1-(1+(Mg/3.57).*exp(-V/V0)).^(-1)))/(V0));
IN_func = @(V, Mg) g*(1+(Mg/3.57).*exp(-V/V0)).^(-1).*(V-EN);

fzero(dIN_func, -40);

V_span = -90:0.001:0;
IN_vals1 = IN_func(V_span, 0);
IN_vals2 = IN_func(V_span, 0.1);
IN_vals3 = IN_func(V_span, 1);


figure (2)
plot(V_span,IN_vals1, 'LineWidth', 1)
x = xlim;
y = ylim;
hold on
%plot(t,u2, 'LineWidth', 1)
%plot(t,u3, 'LineWidth', 1)
%y = ylim; % current y-axis limits
%plot([tpeak1 tpeak1],[y(1) y(2)], 'LineWidth',1) % visual check
plot(V_span,IN_vals2, 'LineWidth', 1)
plot(V_span,IN_vals3, 'LineWidth', 1)
hold off
title('NMDA Receptor Conductance vs. Voltage', 'interpreter','latex','FontSize', 14)
legend('[Mg^{2+}] = 0', '[Mg^{2+}] = 0.1', '[Mg^{2+}] = 1')
xlabel('$V$ (mV)','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$I_N(V)$','interpreter','latex', 'FontWeight','bold','FontSize',16)




%==========================================================================








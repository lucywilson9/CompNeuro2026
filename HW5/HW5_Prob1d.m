%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 5, Problem 1d
% November 2, 2026

%--------------------------------------------------------------------------
%% Wilson-Cowan equations
% Choose wEE = (1+eps)w_EE^Hopf for a small eps>0. For TI = 10, e = 0.005,
% 0.01, 0.02, 0.05, 0.1, 0.2, integrate from a point close to (Ebar, Ibar)
% until the oscillation has settled and report the period, the amplitude
% max E - min Eand the lag of peaks of I behind those of E (in degrees).
% Plot the amplitude against eps on a log-log scale with a fitted slope.
% for TI = 2, 3, 5, 7.5, 10, 15, 20, eps = 0.02, compare the simulated
% frequency with omega/2pi from (b).
% plot E(t), I(t) for TI = 10, eps = 0.1. How do the period and lag change
% as eps grows, and why? Is the Hopf bifurcation supercritical?

% girl what does the w_EE^hopf mean.. I think it's the trace 0 val?

% w_ee^Hopf

eps_vec = [0.005, 0.01, 0.02, 0.05, 0.1, 0.2];
tI_vec = [2, 3, 5, 7.5, 10, 15, 20];
E_amp_vec = zeros(length(eps_vec),1);
period_vec = [];
freq_vec = [];
pred_freq_vec = [];
lag_deg_vec = [];
for i = 1:length(eps_vec) % or length(tI_vec)
eps = eps_vec(i);
%eps = 0.02;
tI = 10;
%tI = tI_vec(i);
% just use an ode solver for a while until period settles
wEE_hopf = 4.167;
wEE = (1+eps)*wEE_hopf;
Ebar = 0.1;
Ibar = 0.1;
wEI = 10;
wIE = 10;
tE = 5;
n = 4;
pI = -1.549; % from part a I think
pE = 1-log(3)/2-0.1*wEE; % from part a I think


tspan = [0 400];
y0 = [Ebar + 0.002 Ibar - 0.002];
[t,y] = ode15s(@(t,y) odefcn(t,y,wEE, wEI, wIE, pE, pI, n, tE, tI), tspan, y0);
E = y(:,1);
I = y(:,2);

% find peaks
E_peak_vec = [];
I_peak_vec = [];
for j = 2:(length(E)-1)
    if E(j+1) < E(j) && E(j) > E(j-1)
        E_peak_vec(end+1) = find(E==E(j));
    end
    if I(j+1) < I(j) && I(j) > I(j-1)
        I_peak_vec(end+1) = find(I==I(j));
    end
end

period_vec(end+1) = t(E_peak_vec(end)) - t(E_peak_vec(end-1));
freq_vec(end+1) = 1/period_vec(end);

pred_freq_vec(end+1) = sqrt((n^2*wEI*wIE*Ebar*(1-Ebar)*Ibar*(1-Ibar))/(tE*tI) - 1/tI^2);
 
E_amp = max(E) - min(E);
E_amp_vec(i) = E_amp;

lag_time = mean(I_peak_vec(1:min(length(I_peak_vec), length(I_peak_vec))) - E_peak_vec(1:min(length(E_peak_vec), length(I_peak_vec))));
lag_deg = mod((lag_time / period_vec(end)) * 360, 360);
if lag_deg > 180, lag_deg = lag_deg - 360; 
end

lag_deg_vec(end+1) = lag_deg;

end

logeps = log10(eps_vec); 
logamp = log10(E_amp_vec);
p = polyfit(logeps, logamp, 1);
slope = p(1); 
intercept = 10^p(2);


xFit = logspace(log10(min(eps_vec)), log10(max(eps_vec)), 100);
yFit = intercept * (xFit.^slope);



% ok, this works - now need the period and lag for this
figure (1)
plot(t,E, 'LineWidth', 1)
x = xlim;
y = ylim;
hold on
plot(t, I, 'LineWidth',1)
hold off
title('$E(t), I(t), \epsilon=0.1$', 'interpreter','latex','FontSize', 14)
legend('E(t)', 'I(t)')
xlabel('time','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('Populations','interpreter','latex', 'FontWeight','bold','FontSize',16)

figure (2)
loglog(eps_vec,E_amp_vec, '*', 'LineWidth', 1)
x = xlim;
y = ylim;
hold on
loglog(xFit, yFit, '-r', 'LineWidth', 2, 'DisplayName', 'Linear Fit');
hold off
title('$Amp_E, \epsilon$ varied', 'interpreter','latex','FontSize', 14)
legend('E(t)', 'I(t)')
xlabel('$\epsilon$','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('Amplitude of E','interpreter','latex', 'FontWeight','bold','FontSize',16)






% =========================================================================
function dydt = odefcn(~,y,wEE, wEI, wIE, pE, pI, n, tE, tI)
dydt = zeros(2,1);
dydt(1) = -y(1)/tE + 1./(tE*(1 + exp(-n.*(wEE*y(1)-wEI*y(2)+pE)))); % E
dydt(2) = -y(2)/tI + 1./(tI*(1 + exp(-n.*(wIE*y(1)+pI)))); % I
end
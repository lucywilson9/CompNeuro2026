%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 3, Problem 4c
% October 4, 2026

%--------------------------------------------------------------------------
%% Measuring Lag with Gap-Junction coupled QIF Neurons
% Then simulate the voltage equations with the infinite threshold and reset replaced by ±L, L = 200, 
% starting neuron 1 at reset and neuron 2 a fraction φ0 = 0.4 of a cycle behind on the uncoupled orbit. 
% Measure the lag once per cycle from spiketimes, φk = (t(2)k − t(1)k)/(t(1)k+1 − t(1)k), and plot it 
% against the prediction for  = 0.02. Fit the decay rate oftan(πφk) for  = 0.005, 0.02, 0.05 and compare with 2

L = 200;
dt = 0.0001;
initial_time_lag = 0.6*pi; % 0.4 of a cycle behind reset = 60% of the way to threshold
eps = 0.02; % tweak as needed

u10 = -L;
u20 = -cot(initial_time_lag); % uncoupled orbit follows v(t) = -cot(t)

nt = round(90/dt) + 1;
t  = linspace(0, 90, nt)';
u1  = zeros(nt,1);
u2  = zeros(nt,1);
u1(1) = u10;
u2(1) = u20;
spikes_1 = [];
spikes_2 = [];

% QIF solve for both neurons:
for j = 1:nt-1
    du1 = dt*(1 + u1(j)^2 + eps*(u2(j)-u1(j)));
    du2 = dt*(1 + u2(j)^2 + eps*(u1(j)-u2(j)));
    u1(j+1) = u1(j) + du1;
    u2(j+1) = u2(j) + du2;
    if u1(j+1) >= L
       u1(j+1) = -L;
       spikes_1(end+1,1) = t(j+1);   %#ok<SAGROW>
    end
    if u2(j+1) >= L
       u2(j+1) = -L;
       spikes_2(end+1,1) = t(j+1);   %#ok<SAGROW>
    end
end

% finding the measured phis
phis = [];
phis_times = [];

for i = 1:(length(spikes_1) - 1)
    % make sure we are correlating the spikes from diff neurons correctly
    index = find(spikes_2 > spikes_1(i) & spikes_2 < spikes_1(i+1), 1, 'first');
    spike_2 = spikes_2(index);

    phik = (spike_2 - spikes_1(i))/(spikes_1(i+1)-spikes_1(i));
    phis = [phis, phik]; %#ok<AGROW>
    phis_times = [phis_times, spike_2]; %#ok<AGROW> % record the time to use w/expected formula, plotting
end

% find expected val for phi by making the implicit hw formula explicit
phi_expec = atan(tan(0.4*pi)*exp(-2*eps*phis_times))/pi;

% "decay rate of tan(pi phi)" if we log both sides we are fitting a linear
% approx.
dat = log(tan(pi * phis)); % put data in an easier (log-linear) form for fitting decay func. 

decay_line = polyfit(phis_times, dat, 1); % fit a linear decay model to the measured phi_ks
fitted_rate = -decay_line(1); % get the decay rate, negate so sign matches sign of 2eps
  
fprintf('eps = %-10.3f\n2eps = %-20.3f\nfitted rate = %-20.3f\n', eps, 2*eps, fitted_rate);

figure (1)
plot(phis_times, phis, '*-', 'LineWidth', 1);
hold on;
plot(phis_times, phi_expec, '-', 'LineWidth', 1);
hold off
legend('Measured \phi_k', 'Theoretical \phi_k')
xlabel('Time (t)','FontWeight','bold','FontSize',14);
ylabel('Phase Lag $\phi$', 'interpreter','latex', 'FontWeight','bold','FontSize',14);
title("Measured vs. Expected Phase Lag, $\epsilon = 0.02$",'interpreter','latex','FontWeight','bold','FontSize', 16)

%==========================================================================



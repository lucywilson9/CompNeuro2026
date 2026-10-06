%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 3, Problem 1c
% October 4, 2026

%--------------------------------------------------------------------------
%% Simulating Mode-locking
% Simulate the model for  = 0.5 and I ∈ [0.9, 3], and plot the firing number (spikes per forcing period,
% measured after a transient) against I, along with the unforced firing number 1/ log[I/(I − 1)] for I > 1
% (zero for I ≤ 1). Mark the tongue boundaries from part (a) and compare with the measured 1:1 plateau.
% Identify at least four other plateaus, and report the phase of the spikes at I = 1.6. Hint: between spikes
% v(t) = G(t) + [v(Tn) − G(Tn)]e−(t−Tn) exactly, so you only need to locate threshold crossings.

eps = 0.5;
I_range = 0.9:0.001:3; % range for constant in forced_lif function
T = 65;
dt = 0.001;

firing_nos = [];
for i=1:length(I_range)
    [spikes, phases] = forced_lif(I_range(i), eps, T, dt);
    if sum(spikes > 0.2*T) >= 2 % "burn-in," don't count spikes as important unless we see at least two to begin with
        firing_nos(end+1) = length(spikes)/T; %#ok<SAGROW
    else 
        firing_nos(end+1) = 0;
    end
end

[spikes_fix, phases_fix] = forced_lif(1.6, eps, T, dt);

firing_no_unforced = @(I) (0).*(I <= 1) + (1./(log(I./(I-1)))).*(I > 1);

firing_no_unforceds = firing_no_unforced(I_range);
left_bdry = 1/(1-exp(-1)) - eps/(sqrt(4*pi^2+1));
right_bdry = 1/(1-exp(-1)) + eps/(sqrt(4*pi^2+1));

figure(1)
plot(I_range, firing_nos)
y = ylim;
x = xlim;
hold on
plot(I_range, firing_no_unforceds)
plot([left_bdry left_bdry], [y(1) y(2)], 'r--', 'LineWidth',1)
plot([right_bdry right_bdry], [y(1) y(2)], 'b--', 'LineWidth',1)
hold off
xlabel("Input $I$", 'interpreter','latex','FontWeight','bold','FontSize',14)
ylabel("No. Spikes/Forcing Period", 'FontWeight','bold','FontSize',14)
legend("Measured Firing Numbers", "Expected Firing Numbers", 'Left tongue bdry, $\sim 1.5034', 'Right tongue bdry, $\sim 1.6606')
title("Firing Rates vs. Periodic Forcing Input currents, $T=65$",'interpreter','latex','FontWeight','bold','FontSize', 16)






%==========================================================================
function [spikes, phases] = forced_lif(I, eps, T, dt) % code from spike_map.m
nt = round(T/dt) + 1;
t  = linspace(0, T, nt);
u  = 0; spikes = [];
for j = 1:nt-1
    u = u + dt*(I + eps*sin(2*pi*t(j)) - u);
    if u > 1
        u = 0; spikes(end+1,1) = t(j+1);   %#ok<AGROW>
    end
end
phases = mod(spikes, 1);
end
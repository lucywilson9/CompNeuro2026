%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 6, Problem 1c
% November 6, 2026

%--------------------------------------------------------------------------
%% Foraging over patches
% (i) Simulate 4 × 10^5 patch visits at T = 2, T^∗ and 10 and compare the 
% measured intake with R(T). (ii) Plot R(T) together with the line Rknown 
% and the curve λ˜(T), marking T^∗ and R^∗.

% params:
p = 0.5;
ph = 0.5;
tau = 5;
num_patches = 4*10^5;

R_T = @(T) (ph*(1-exp(-p.*T)))./(tau+(1-ph).*T+(ph*(1-exp(-p.*T)))/p);
lambda_T = @(T) (p.*ph.*exp(-p.*T))./(ph.*exp(-p.*T)+(1-ph));
R_known = ph/(tau+ph/p);

Topt = fzero(@(T) R_T(T) - lambda_T(T), 4.0); % optimal time
Ropt= R_T(Topt); % optimal rate
T_vals = [2,Topt,10];

% patch visits sim
rng(9); % set random seed

fprintf('%-10s %-15s %-15s %-15s\n', 'T', 'Measured R', 'Analytical R');
for i = 1:length(T_vals)
    T = T_vals(i);
    high_vec = rand(num_patches,1) < ph; % vector assigning each patch to be high or low with prob ph
    dinner_time_vec = exprnd(1/p,num_patches,1); % vector of exp drawn discovery times
    
    eat = high_vec&(dinner_time_vec<=T); % animal gets the chunk if they visit a high patch and find it before wait time T
    
    time_at_patch = zeros(num_patches, 1); % patch time: min(exp draw, T)
    time_at_patch(high_vec) = min(dinner_time_vec(high_vec),T);
    time_at_patch(~high_vec) = T; % if it was a low patch, set the time at patch to t
    
    R_meas = sum(eat)/sum(time_at_patch+tau); % overall rate (add travel time to each)
    R_analyt = R_T(T);
    
    fprintf('\n%-10.4f %-15.5f %-15.5f', T, R_meas, R_analyt);
end

% plotting stuff
T_grid = linspace(0, 10, 500);

figure('Color', 'w')
plot(T_grid, R_T(T_grid), 'b', 'LineWidth', 2)
hold on
plot(T_grid, lambda_T(T_grid), 'r', 'LineWidth', 2)
yline(R_known, 'k', 'LineWidth', 2)
plot(Topt, Ropt, 'co', 'MarkerFaceColor', 'c', 'MarkerSize', 8)
hold off
legend('R(T)', '\lambda(T)', 'R_{known}','(T^*, R^*)')
xlabel('Wait Time T')
ylabel('Rate R(T)')
title('Overall vs. Current Food Rate, Varying Wait Time')



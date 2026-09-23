%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 2, Problem 4c
% September 22, 2026

%--------------------------------------------------------------------------
%% Simulating Shot noise
% fix <V> = 10mV = wvtau, τ = 20 ms fixed while taking ντ = 2, 20, 200 with w set accordingly
% Report the measured mean, variance, and skewness against theory in a table
% for vt = 2, 200:overlay the voltage histogram with the matched Gaussian

tau = 0.020; %ms
vtau = 200;
v = vtau/tau;
w = 10/vtau;

T = 500; % arb for now, ms
dt = 0.001;
times = 0:dt:T;


volt_curve = volt_func(times, tau, w, dt, v); % skeptical this will work right away
volt_curve = volt_curve(1/10*length(curve-1):end);
report_mean = mean(volt_curve);
%report_std = std(volt_curve);
report_var = var(volt_curve);
report_skew = skewness(volt_curve);

theor_skew = 2*sqrt(2)/(3*sqrt(vtau));
theor_var = v*w^2*tau/2;

fprintf('meas. mean: %.13f\nmeas. var: %.13f\nmeas. skew: %.13f\n', report_mean, report_var, report_skew)
fprintf('theor. mean: 10mV\ntheor. var: %.13f\ntheor. skew: %.13f\n', theor_var, theor_skew)

% for plotting the gaussian:
mu = report_mean;
std = sqrt(report_var);
x = (mu-4*std):0.001:(mu+4*std);
gauss = normpdf(x, mu, std);

figure (3)
histogram(volt_curve,'Normalization','pdf')
hold on
plot(x, gauss, 'LineWidth',1)
hold off
xlabel("Voltage (mV)", 'FontWeight','bold','FontSize',14)
ylabel("pdf", 'FontWeight','bold','FontSize',14)
legend("V(t)", "Norm(<V>, sqrt(Var[V]))")
title("Voltage Histogram: $\nu\tau=200$",'interpreter','latex','FontWeight','bold','FontSize', 16)


%==========================================================================
function volt_vec = volt_func(t, tau, w, dt, v)
volt_vec =[];
volt_vec(1) = 0;
for j = 2:length(t)
        volt_vec(end+1) = volt_vec(end)*exp(-dt/tau) + w*poissrnd(v*dt);

end
end







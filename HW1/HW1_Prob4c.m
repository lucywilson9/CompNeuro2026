%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 1, Problem 4c
% September 9, 2026

%--------------------------------------------------------------------------
%% Write code to integrate the LIF equation with this input and verify the formula from (b).  
% define some params
R = 1;
tm = 10;
ts = 2;
urest = 0;
uth = 1;
r = ts/tm;

I0 = (uth/R)*(r)^(1/(r-1));
I0b = 0.75*I0;
I0a = 1.25*I0;

I0_crit = @(uth, R, ts, tm) (uth/R)*(ts./tm).^(1./((ts./tm)-1));

ts_span = 0.5:0.001:100;
take_out= 10;
ts_span_no_div = setdiff(ts_span, take_out);
I0_crit_vals = I0_crit(uth, R, ts_span_no_div, tm);

Ivec = [I0 I0b I0a]; 

T = 50; % arb endpoint
dt = 0.1;

tspan = 0:0.001:T;
u0 = urest;

% note that I am aware that there are more efficient ways of executing
% these
uprime1 = @(t,u) (1/tm)*(-u + urest + R*I0*exp(-t./ts));
uprime2 = @(t,u) (1/tm)*(-u + urest + R*I0b*exp(-t./ts));
uprime3 = @(t,u) (1/tm)*(-u + urest + R*I0a*exp(-t./ts));

[t,u1] = ode45(uprime1, tspan, u0);
[t,u2] = ode45(uprime2, tspan, u0);
[t,u3] = ode45(uprime3, tspan, u0);

tvals = 0:0.1:T;
tpeak = ((ts*tm)/(ts-tm))*log(ts/tm);
umax1 = R*I0*(r)^(-1/(r-1));
umax2 = R*I0b*(r)^(-1/(r-1));
umax3 = R*I0a*(r)^(-1/(r-1));

tpeak1 = t(find(u1==max(u1)));
tpeak2 = t(find(u2==max(u2)));
tpeak3 = t(find(u3==max(u3)));


fprintf('\nI0 crit umax: %.4f\nI0 small umax: %.4f\nI0 big umax: %.4f\n', max(u1), max(u2), max(u3));
fprintf('\nI0 crit: %.4f\nI0 small: %.4f\nI0 big: %.4f\n', I0, I0b, I0a);
fprintf('\nI0 crit t_peak: %.4f\nI0 small t_peak: %.4f\nI0 big t_peak: %.4f\n', tpeak1, tpeak2, tpeak3);

figure (1)
plot(t,u1, 'LineWidth', 1)
x = xlim;
y = ylim;
hold on
plot(t,u2, 'LineWidth', 1)
plot(t,u3, 'LineWidth', 1)
y = ylim; % current y-axis limits
%plot([tpeak1 tpeak1],[y(1) y(2)], 'LineWidth',1) % visual check
plot([x(1) x(2)], [uth uth], 'k--', 'LineWidth',1)
hold off
title('LIF Model', 'FontSize', 14)
legend('I_0 = I_0^{crit}', 'I_0 < I_0^{crit}', 'I_0 > I_0^{crit}', 'u_{th}')
xlabel('$t$','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$u(t)$','interpreter','latex', 'FontWeight','bold','FontSize',16)


figure (2)
loglog(ts_span_no_div, I0_crit_vals, 'LineWidth',1)
title('I0 vs. tau_s', 'FontSize', 14)
title('$I_0^{crit}$ vs. $\tau_s$','interpreter','latex', 'FontWeight','bold','FontSize',16)
xlabel('$\tau_s$','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$I_0$','interpreter','latex', 'FontWeight','bold','FontSize',16)



%==========================================================================








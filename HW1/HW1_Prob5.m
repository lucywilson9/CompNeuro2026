%% APPM 5370 Computational Neuroscience
% Lucy Wilson
% Homework 1, Problem 5
% September 9, 2026

%--------------------------------------------------------------------------
%% Exploring the rheobase in the HH Model

Id = 7;
tmax = 300;
gK = 36;
gNa = 120;
skipfrac = 0.4;

% need len

% first: finding the rheobase with bisection
% Id is the constant current, so rheobase is the smallest constant current
% that gives sustained repetitive firing.
Id_range = [5 7]; % ad hoc/comments informed.

root = bisection(5, 7, 10^(-16), 0);
fprintf('rheobase val: %.13f\n', root)
Id_rheo = root;
Id_rheo_b = root - 10^(-2);
Id_rheo_a = root + 10^(-2);

root = bisection(1, 4, 10^(-16), 1);
fprintf('rheobase val: %.13f\n', root)

rheo_seq = Id_rheo + 10^(-12):0.1:2*Id_rheo;
rate_vec = [];

rheo_seq2 = root + 10^(-12):0.1:2*root;
rate_vec2 = [];

for i = 1:length(rheo_seq)
    [t, V, spk, rate] = hh_sim(rheo_seq(i), 300, 36, 120, 0.4);
    rate_vec(end+1) = rate;
end

for i = 1:length(rheo_seq2)
    [t, V, spk, rate2] = hh_sim(rheo_seq2(i), 300, 30, 120, 0.4);
    rate_vec2(end+1) = rate2;
end

% check: currently reporting rheobase of 6.257

[t, V, spk, rate] = hh_sim(6.257, 300, 36, 120, 0.4); % does yield length(spk) = 2

% now need to plot stuff against each other.
[t1, V1, spk, rate] = hh_sim(Id_rheo, 300, 36, 120, 0.4);
[t2, V2, spk, rate] = hh_sim(Id_rheo_b, 300, 36, 120, 0.4);
[t3, V3, spk, rate] = hh_sim(Id_rheo_a, 300, 36, 120, 0.4);

% restrict to after the 'skipfrac' period
t1_to_plot = t1(find(t1==120):end);
V1_to_plot = V1(find(t1==120):end);
t2_to_plot = t2(find(t2==120):end);
V2_to_plot = V2(find(t2==120):end);
t3_to_plot = t3(find(t3==120):end);
V3_to_plot = V3(find(t3==120):end);

% figure(3)
% plot(t1_to_plot, V1_to_plot, 'LineWidth',1)
% hold on
% plot(t2_to_plot, V2_to_plot, 'LineWidth',1)
% plot(t3_to_plot, V3_to_plot, 'LineWidth',1)
% hold off
% legend('Id = rheobase', 'Id < rheobase', 'Id > rheobase')
% title('Hodgkin-Huxley simulations, Id near rheobase')
% xlabel('$t$','interpreter','latex', 'FontWeight','bold','FontSize',16)
% ylabel('$V(t)$','interpreter','latex', 'FontWeight','bold','FontSize',16)

figure (4)
plot(rheo_seq, rate_vec,'r', 'LineWidth',1)
x = xlim;
hold on
plot(rheo_seq2, rate_vec2, 'b','LineWidth',1)
y = ylim;
plot([Id_rheo Id_rheo], [y(1) y(2)], 'r--', 'LineWidth',1)
plot([root root], [y(1) y(2)], 'b--', 'LineWidth',1)
hold off
legend('gK=36', 'gK=30', 'rheobase 1 ~ 6.257', 'rheobase 2 ~ 2.773')
title('$f-I$ plot, rheobase to double rheobase','Interpreter','latex','FontWeight', 'bold', 'FontSize',16)
xlabel('$I$','interpreter','latex', 'FontWeight','bold','FontSize',16)
ylabel('$f$','interpreter','latex', 'FontWeight','bold','FontSize',16)
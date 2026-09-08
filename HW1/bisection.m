%% APPM 5470 Computational Neuroscience
% Lucy Wilson
% Homework 1, Problem 5
% September 9, 2026

%--------------------------------------------------------------------------
%% Bisection method
% Note: following pseudocode from Sauer 2005 text

function root = bisection(a,b,tol, ind)

if ind == 1
    gK = 30;
else
    gK = 36;
end

[t, V, spka, rate] = hh_sim(a, 300, gK, 120, 0.4);
[t, V, spkb, rate] = hh_sim(b, 300, gK, 120, 0.4);
fa = length(spka) - 2;
fb = length(spkb) - 2;
if sign(fa)*sign(fb) >= 0 
    return % if both approximations for y(1) are of the same sign, no 
    % root in the interval
end

%fa = f(a);
%fb = f(b);

while (b-a)/2 > tol % while width of interval is greater than tolerance
%while (b-a) > tol
    c = (a+b)/2; % define midpt
    [t, V, spkc, rate] = hh_sim(c, 300, gK, 120, 0.4);
    fc = length(spkc) - 2;
    %fc = f(c); 
    if fc ==0
        break
    end
    if sign(fc)*sign(fa) < 0 
        b = c; % save appropriate endpoint of region with root
        fb = fc;
    else
        a = c;
        fa = fc;
    end
end
root = (a+b)/2;

end
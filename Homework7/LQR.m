%% Reset
clear variables
close all
clc

%% Part 1
R_values = logspace(-3,log10(2)/log10(10),100);
CLeigs = zeros(2,numel(R_values));
max_rise = 0;
min_rise = Inf;
max_overshoot = 0;
min_overshoot = Inf;
for i = 1:numel(R_values)
tf_num = 6543;
tf_denom = [1 75.74 3137];
tf_sys = tf(tf_num,tf_denom);
ss_canon = canon(tf_sys,'companion');
Actrb = ss_canon.A';
Bctrb = ss_canon.B;
Cctrb = ss_canon.C;
Dctrb = ss_canon.D;

%% Part 2
R = R_values(i);
% R = 1;
Q = Cctrb'*Cctrb;

%% Part 3
[K,S,E] = lqr(ss_canon.A,ss_canon.B,Q,R);
CLeigs(:,i) = E;

%% Part 4
obsv_lam = [3*real(E(1)) 5*real(E(1))];

%% Part 5
mu = zeros(3,2);
v = zeros(2,2);
theta = zeros(1,2);
%null space is mu
for j = 1:numel(obsv_lam)
    S = [obsv_lam(j)*eye(2)-Actrb',Cctrb'];
    [U,S,V] = svd(S);
    mu(:,j) = V(:,3);
    v(:,j) = mu(1:2,j);
    theta(:,j) = mu(3,j);
end
H = (-theta/v)';

%% Part 6
compensator_ss = ss(Actrb-Bctrb*K-H*Cctrb,H,K,0);
compensator_tf = tf(compensator_ss);
a = step(compensator_tf);
b = stepinfo(a);
if b.RiseTime > max_rise
    max_rise = b.RiseTime;
end
if b.RiseTime < min_rise
    min_rise = b.RiseTime;
end
if b.Overshoot > max_overshoot
    max_overshoot = b.Overshoot;
end
if b.Overshoot < min_overshoot
    min_overshoot = b.Overshoot;
end
if b.RiseTime <= 0.01 && b.Overshoot <= 5
    figure(1)
    plot(1:numel(a),a)
    break
end
end
plot(real(CLeigs),imag(CLeigs), 'kx')
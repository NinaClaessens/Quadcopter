close all;
clear

load('references_16.mat');

m = 0.5;
L = 0.25;
k = 3*10^-6;
b = 10^-7;
g = 9.81;
kd = 0.25;
Ixx = 5*10^-3;
Iyy = 5*10^-3;
Izz = 1*10^-2;
cm = 10^4;
uhover = (g*m)/(4*cm*k);
u1 = uhover;
u2 = uhover;
u3 = uhover;
u4 = uhover;

% construct Jacobian = A
A= zeros(12);
A(1,4) = 1;
A(2,5) = 1;
A(3,6) = 1;
A(4,4) = -kd/m;
A(4,8) = k*cm/m*(u1+u2+u3+u4);
A(5,5) = -kd/m;
A(5,7) = -k*cm/m*(u1+u2+u3+u4);
A(6,6) = -kd/m;
A(7,10) = 1;
A(8,11) = 1;
A(9,12) = 1;

% Jacobian for inputs = B
B = zeros(12,4);
B(6,1:4) = k*cm/m;
B(10,1) = L*k*cm/Ixx;
B(10,3) = -L*k*cm/Ixx;
B(11,2) = L*k*cm/Iyy;
B(11,4) = -L*k*cm/Iyy;
B(12,1) = b*cm/Izz;
B(12,2) = - b*cm/Izz;
B(12,3) = b*cm/Izz;
B(12,4) = - b*cm/Izz;

% construct C
C = eye(12);

D = zeros(12,4);

sys = ss(A,B,C,D);

% discretization
Ts = 0.05;
[Ad, Bd, Cd, Dd] = bilinear(A,B,C,D,1/Ts);
sysd = ss(Ad,Bd,Cd,Dd,Ts);

% LQR controller
Qposition = [1000 1000 1000];
Qvelocity = [0.1 0.1 0.1];
Qangle    = [0.1 0.1 0.1];
Qangveloc = [0.1 0.1 0.1];
Q = diag([Qposition Qvelocity Qangle Qangveloc]);
R = eye(4);
K = dlqr(Ad,Bd,Q,R);

% N = ones(12,4);
% K = dlqr(Ad,Bd,Q,R);

K = lqr(sys,Q,R);

N =pinv([Ad-eye(12), Bd; Cd, Dd])*[zeros(12,12); eye(12)];
Nx = N(1:12,:);
Nu = N(13:end,:);
% probleem: unobservable mode on unit circle (denk ik) ... 
% https://nl.mathworks.com/help/control/ref/dlqr.html
% bij limitations
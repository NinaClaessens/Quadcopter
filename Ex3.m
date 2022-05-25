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
sysd = c2d(ss(A,B,C,D),Ts,'tustin');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;
% [Ad, Bd, Cd, Dd] = bilinear(A,B,C,D,1/Ts); % zelfde als tustin



% LQR controller
Qposition = [100 100 100];
Qvelocity = [1 1 1];
Qangle    = [1000 1000 1000];
Qangveloc = [1 1 1];
Q = diag([Qposition Qvelocity Qangle Qangveloc]);
R = eye(4)*0.01;
% K = dlqr(eye(12)+A*Ts,B*Ts,Q,R); euler
K = dlqr(Ad,Bd,Q,R);

sysd = ss(Ad,Bd,Cd,Dd,Ts);

% N =pinv([A*Ts, B*Ts; C, D])*[zeros(12,12); eye(12)]; %euler
N =pinv([Ad-eye(12), Bd; Cd, Dd])*[zeros(12,12); eye(12)];
Nx = N(1:12,:);
Nu = N(13:end,:);

% Nbar = rscale(sysd,K);

%% integral control
% augmented system
Aa = [eye(size(Cd,1)) Cd; zeros(size(Ad,1),size(Cd,1)) Ad];
Ba = [Dd; Bd];
Co = ctrb(Aa,Ba);
disp('nb of uncontrollable states is');
disp(length(Aa) - rank(Co))

Q_int = [zeros(size(Q)) Q; Q zeros(size(Q))];
K_int = lqi(sysd,Q_int,R); % geeft zelfde error
% K_int = dlqr([eye(12), Cd; zeros(12), Ad],[Dd;Bd], eye(24),eye(4));
% K1 = K_int(1:12,:);
% K0 = K_int(13:24,:);
% probleem: unobservable mode on unit circle (denk ik) ... 
% https://nl.mathworks.com/help/control/ref/dlqr.html
% bij limitations
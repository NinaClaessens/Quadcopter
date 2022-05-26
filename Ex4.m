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
C = zeros(3,12);
C(1:3,1:3) = eye(3);
D = zeros(3,4);

% continuous system
sys = ss(A,B,C,D);


Ts = 0.05;
% Ad = A*Ts+eye(12); %euler
% Bd = B*Ts;
% Cd = C;
% Dd = D;
[Ad, Bd, Cd, Dd] = bilinear(A,B,C,D,1/Ts); % bilinear
sysd = ss(Ad,Bd,Cd,Dd,Ts);

Aa = [eye(3), Cd; zeros(12,3), Ad];
Ba = [Dd;Bd];
Ca = [zeros(3) Cd];

rank(ctrb(Aa, Ba)) %controllable!
Qerr = [10 10 10];
Qposition = [700 700 4000]; 
Qvelocity = [0.1 0.1 0.1];
Qangle    = [2500 3500 4000];
Qangveloc = [1 1 1];
AugQ = diag([Qerr Qposition Qvelocity Qangle Qangveloc]);
AugR = eye(4)*3.4/200;
K_int = dlqr(Aa,Ba,AugQ,AugR); 
K1 = K_int(:,1:3);
K0 = K_int(:,4:15);

%% Kalmann
sys_k = ss(Ad,[Bd eye(12)],Cd,[Dd zeros(3,12)],Ts); % met noise inputs
rank(obsv(Ad,Cd)) %niet observable
Q_k = zeros(12);
Q_k(1:3,1:3) = eye(3)*2.5*10^(-5);
Q_k(7:9,7:9) = eye(3)*7.57*10^(-5);
R_k= eye(3)*2.5*10^(-5); %=random nu, moet denk ik measurement noise op y zijn
kalmf = kalman(sys_k,Q_k,R_k)




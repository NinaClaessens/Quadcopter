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


Ts = 0.05;

% bilinear
sysd = c2d(sys,Ts,'tustin');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;


%% LQR controller
Qposition = [70 70 9000]; 
Qvelocity = [0.1 0.1 0.1];
Qangle    = [850 850 50];
Qangveloc = [1 1 1];
Q = diag([Qposition Qvelocity Qangle Qangveloc]);
R = eye(4)*3.4;

K = dlqr(Ad,Bd,Q,R);

%pole plot
Acl = Ad-Bd*K;
Bcl = Bd*K;
Ccl = Cd-Dd*K;
Dcl = Dd*K;

syscl = ss(Acl,Bcl,Ccl,Dcl);
theta = 0:0.01:2*pi;
figure
pzplot(syscl)
hold on
plot(cos(theta),sin(theta))

N =pinv([Ad-eye(12), Bd; Cd, Dd])*[zeros(12,12); eye(12)];
Nx = N(1:12,:);
Nu = N(13:end,:);

%% Integral controller
% construct C
C = zeros(3,12);
C(1:3,1:3) = eye(3);
D = zeros(3,4);
sysd = c2d(ss(A,B,C,D),Ts,'tustin');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

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
AugR = R/200;
K_int = dlqr(Aa,Ba,AugQ,AugR); 
K1 = K_int(:,1:3);
K0 = K_int(:,4:15);

%generate pole plot
Acl = Aa-Ba*K_int;
Bcl = Ba*K_int;
Ccl = Ca-Dd*K_int;
Dcl = Dd*K_int;
syscl = ss(Acl,Bcl,Ccl,Dcl);
theta = 0:0.01:2*pi;
figure
pzplot(syscl)
hold on
plot(cos(theta),sin(theta))
clear
close all

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
u1 = (g*m)/(4*cm*k);
u2 = (g*m)/(4*cm*k);
u3= (g*m)/(4*cm*k);
u4 = (g*m)/(4*cm*k);

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
C = zeros(6,12);
C(1,1) = 1;
C(2,2) = 1;
C(3,3) = 1;
C(4,7) = 1;
C(5,8) = 1;
C(6,9) = 1;

% construct D
D = zeros(6,4);

% sampling time
Ts = 0.05;

%% zero order hold
% Ad = exp(A*Ts);
% Bd = (Ad-eye(size(A)))*A\B;       % werkt niet want A singulier
% Cd = C;
% Dd = D;

%% Euler's method
% Ad = eye(size(A)) + Ts*A;
% Bd = Ts*B;
% Cd = C;
% Dd = D;

%% Bilinear transformation
sysd = c2d(ss(A,B,C,D),Ts,'tustin');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

%% Analysis

% stability
[U,V] = eig(Ad);
disp('eigenvalues are')
disp(diag(V))     % not strictly within unit circle 
                  % => no internal stability

% poles
disp('poles are')
disp(pole(sysd))   % same as eigenvalues
                   % => minimal system
                   % => not in/out stable
disp('zeros are')
disp(tzero(sysd))

% controllabililty
Co = ctrb(Ad,Bd);
disp('nb of uncontrollable states is');
disp(length(Ad) - rank(Co))     % everything controllable -> stabilizable

% observability
Ob = obsv(Ad,Cd);
disp('nb of unobservable states is');
disp(length(A)-rank(Ob))        % everything observable -> detectable


%% simulation data
t = 0:0.05:5;
dr = zeros(length(t),5);
dr(:,1) = t;
dr(20:end,2:5) = 5;
dr(1:20,2:5) = 0;

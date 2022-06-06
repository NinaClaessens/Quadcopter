% close all;
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

%% Choose poles for LQR 
poles = zeros(12,1);
% select dominant poles first continuous and then transform to discrete.
damping = 0.86; % generally good choice
alpha = -1.05; % close to imaginary axis
ratio_ab = -damping/sqrt(1-damping^2);
beta = alpha/ratio_ab;
poles(1) = exp((alpha + 1i*beta)*Ts);
poles(2) = exp((alpha - 1i*beta)*Ts);


damping = 0.96; % higher for non dominant poles
ratio_ab = -damping/sqrt(1-damping^2);
alpha = alpha*3:(alpha*2-alpha*3)/4:alpha*2;
beta = alpha/ratio_ab;

for i = 1:5
    poles(2*i+1) =  exp((alpha(i) + 1i*beta(i))*Ts);
    poles(2*i+2) =  exp((alpha(i) - 1i*beta(i))*Ts);
end
theta = 0:0.01:2*pi;
figure(1)
plot(real(poles),imag(poles),'o')
hold on
plot(cos(theta),sin(theta))

K_pole = place(Ad,Bd,poles);

% pole plot
Acl = Ad-Bd*K_pole;
Bcl = Bd*K_pole;
Ccl = Cd-Dd*K_pole;
Dcl = Dd*K_pole;

syscl = ss(Acl,Bcl,Ccl,Dcl);
theta = 0:0.01:2*pi;
figure(2)
pzplot(syscl)
hold on
plot(cos(theta),sin(theta))
xlim([-1,1])
ylim([-1,1])

N =pinv([Ad-eye(12), Bd; Cd, Dd])*[zeros(12,12); eye(12)];
Nx = N(1:12,:);
Nu = N(13:end,:);

%% Estimator design

C = zeros(6,12);
C(1:3,1:3) = eye(3);
C(4:6,7:9) = eye(3);
D = zeros(6,4);
sysd = c2d(ss(A,B,C,D),Ts,'tustin');
Ad = sysd.A;
Bd = sysd.B;
Cd = sysd.C;
Dd = sysd.D;

poles = zeros(12,1);
% select dominant poles first continuous and then transform to discrete.
damping = 0.5; % generally good choice
alpha = -0.5; % close to imaginary axis
ratio_ab = -damping/sqrt(1-damping^2);
beta = alpha/ratio_ab;
poles(1) = exp((alpha + 1i*beta)*Ts);
poles(2) = exp((alpha - 1i*beta)*Ts);


damping = 0.8; % higher for non dominant poles
ratio_ab = -damping/sqrt(1-damping^2);
alpha = alpha*3:(alpha*2-alpha*3)/4:alpha*2;
beta = alpha/ratio_ab;

for i = 1:5
    poles(2*i+1) =  exp((alpha(i) + 1i*beta(i))*Ts);
    poles(2*i+2) =  exp((alpha(i) - 1i*beta(i))*Ts);
end
theta = 0:0.01:2*pi;
figure
plot(real(poles),imag(poles),'o')
hold on
plot(cos(theta),sin(theta))

L_pole = place(Ad',Cd',poles)';

sim('pole_placement_controller.mdl')
generate_report(1)

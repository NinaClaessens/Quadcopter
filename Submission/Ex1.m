clear
close all

load('references_16.mat')

% find fixed point
syms k cm m g L Ixx Iyy Izz b u1 u2 u3 u4
eq1 = k*cm/m*(u1 + u2 + u3 + u4) - g == 0;
eq2 = L*k*cm/Ixx*(u1 - u3) == 0;
eq3 = L*k*cm/Iyy*(u2 - u4) == 0;
eq4 = b*cm/Izz*(u1-u2 + u3-u4) == 0;
solve([eq1,eq2,eq3,eq4], [u1,u2,u3,u4])


% define parameters
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
u_hover = (g*m)/(4*cm*k);
u1 = u_hover;
u2 = u_hover;
u3= u_hover;
u4 = u_hover;


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


% short simulation to see difference between real and
% linearized version
t = 0:0.05:5;
dr = zeros(length(t),5);
dr(:,1) = t;
dr(20:end,2:5) = 5;
dr(1:20,2:5) = 0;

r = zeros(length(t),5);
r(:,1) = t;
r(20:end,2:5) = 5+u_hover;   % add hovering voltage
r(1:20,2:5) = 0+u_hover;

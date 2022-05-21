clear all
close all

load('references_16.mat')

%find fixed point
syms k cm m g L Ixx Iyy Izz b u1 u2 u3 u4
eq1 = k*cm/m*(u1 + u2 + u3 + u4) - g == 0;
eq2 = L*k*cm/Ixx*(u1 - u3) == 0;
eq3 = L*k*cm/Iyy*(u2 - u4) == 0;
eq4 = b*cm/Izz*(u1-u2 + u3-u4) == 0;
solve([eq1,eq2,eq3,eq4], [u1,u2,u3,u4])

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

%construct Jacobian
A= zeros(12);
A(1,4) = 1;

A(2,5) = 1;

A(3,6) = 1;

A(4,4) = -kd/m;
A(4,8) = k*cm/m*(u1+u2+u3+u4);

A(5,5) = -kd/m;
A(5,7) = -k*cm/m*(u1+u2+u3+u4);

A(6,6) = -kd/m;

A(7,9) = 1;

A(8,11) = 1;

A(9,12) = 0;

%Jacobian for inputs
B = zeros(12,4);

B(6,1:4) = k*cm/m;

B(10,1) = L*k*cm/Ixx;
B(10,3) = -L*k*cm/Ixx;

B(11,2) = L*k*cm/Iyy;
B(11,4) = -L*k*cm/Iyy;

B(11,1) = b*cm/Izz;
B(11,2) = - b*cm/Izz;
B(11,3) = b*cm/Izz;
B(11,4) = - b*cm/Izz;

%construct c
C = zeros(6,12);
C(1,1) = 1;
C(2,2) = 1;
C(3,3) = 1;
C(4,7) = 1;
C(5,8) = 1;
C(6,9) = 1;

D = zeros(6,4);

t = 0:0.05:10;
r = zeros(length(t),5);
r(:,1) = t;
r(100:end,2:5) = 100 - u1; %substract reference values
r(1:100,2:5) = -u1;

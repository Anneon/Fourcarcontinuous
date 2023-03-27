close all; clc
clear

A = [0 1 0 0 0 0 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 1 0 0;
     0 0 0 0 0 0 0 0;
     0 0 0 0 0 0 0 1;
     0 0 0 0 0 0 0 0];
alpha = 2;
S = 1;
B = [0 0 0 0;
     1 0 0 0;
     0 0 0 0;
     0 1 0 0;
     0 0 0 0;
     0 0 1 0;
     0 0 0 0;
     0 0 0 1];

C = [1 -1 0 0 0 0 0 0;
    0 0 1 -1 0 0 0 0;
    0 0 0 0 1 -1 0 0;
    0 0 0 0 0 0 1 -1];
D = zeros(size(C,1),size(B,2));
states = {'x1', 'h1','x2','h2','x3','h3','x4','h4'};
%inputs = {'ux1', 'uh1','ux2','uh2','ux3','uh3','ux4','uh4'};
inputs = {'u1','u2','u3','u4'};
outputs = {'p1','p2','p3','p4'};
sysc = ss(A,B,C,D);
%step(sysc);

% Define desired eigenvalues for the closed-loop system
des_eigs = [-1 -2 -3 -4 -5 -6 -7 -8];

% Design a state feedback controller
K = place(A,B,des_eigs);

% Define the simulation time and initial state
t = 0:0.01:10;
x0 = [ones(8,1); zeros(8,1)];

% Simulate the closed-loop system
sys = ss(A-B*K, B, C, 0);

%t = 0:0.001:0.05;
sys_cl = ss(A-B*K,B,C,D);
step(sys_cl,t)

%Lqr
Q = 0.05*eye(8);
%R = eye(4);
R=diag([1 10 100 1000]);
U2 = lqr(A-B*K,B,Q,R);

%Observer design

%parameters G and H
G = eye(4); %because 4 states so 4 disturbances (p = n)
H = zeros(4,4); %4 outputs so 0 matrix of 4 x 4

%covariance matrices, process Q, measurement R
Qcov = diag(0.00001*ones(1,4));%Q is 4x4
Rcov = diag(1000*ones(1,4));%R is 4x4
%Bnew = [Bol*G]
sys_kf = ss(A,B*G,C,D*H);
%Obtain Kalman gain L and P,assume W and v are uncorrelated
[kest,L,P] = kalman(sys_kf,Qcov,Rcov,0);
%check the value of L that matlab returns
%Compare L_bar to L(Should be equal)
L_gain = (A*P*C')/(C*P*C'+Rcov);
Error = norm(abs(L_gain-L));
%Assess the stability
Acb = A-L*C;
eig(Acb);

%% Parameters for Optimization
clear all ; clc ; close all ; 
Ts = 1; % Sampling Time  
Ho = 10; % Optimization Horizon (in number of samples)
%% Simulation Time Definition
Tmax = 1e3;
t = 0:Ts:Tmax;
%% System to Control
s = tf('s');
load("tf_delays.mat")
load("realSystem.mat")

sys_d=c2d(Pspm,Ts);
d2d(sys_id_d, Ts);
sys_ss = ss(sys_id_d);

A_Real = sys_ss.A;
B_Real = sys_ss.B;
C_Real = sys_ss.C;
D_Real = sys_ss.D;

%sys_d=c2d(Pspm,Ts)
d2d(sys_id_d, Ts);
sys_ss = ss(sys_id_d);

A = sys_ss.A;
B = sys_ss.B;
C = sys_ss.C;
D = sys_ss.D;

Ly = size(C,1); % Number of lines in the output vector
Lu = size(B,2); % Number of lines in the input vector
Lz = size(A,2); % Number of lines in the state vector

%% Generation of Optimization Matricies 

%INPUTS: 
%feed rate
%mill feed water flow rate
%sump dilution water flow rate
%pump speed

%OUTPUTS:  
%[product particle size  ,  mill solids concentration  ,circulating load,  sump level]

% Constraints on Input/Output
umax = [80 15 45 50]';
umin = [0 2 5 35]';
  
ymax = [100 100 350 2.4]';
ymin = [0       0      0       0]';

% Yref Definition
yref = ymax/2 * ones(1,length(t));
yref_opt = reshape(yref',Ly*length(yref),1);
y0 =ymax/2;

% Reconstruction of yref to make the reference 
% signal compatible with the optimization algorithm

% Weight Matricies
Q = diag([10000 100 1 1]); % Weight on yref - y
R = eye(Lu); % Weight on u(i+1) - u(i)

[M,P,V,H,f,T] = ConstantMatrices(A,B,C,D,Q,R,Ho-1,Ly,Lu);
save ConstantMatrices M P V H f T
%% MPC simulation 
load ConstantMatrices.mat

u = zeros(Lu,length(yref)); % Matrice that allows us to record the consecutive values of the input 
u(1:end ,1 ) = (umax+umin)/2 ;  
% umax = [80 15 45  50]';

y = zeros(Ly,length(yref)); % Mat2rice that allows us to record the consecutive values of the output
y(1:end ,1 ) = y0;

z_Real(: ,1) = pinv(C_Real)*y(1:end ,1) ;   %inv((eye(size(A))-A))*(B)*u(1:end ,1) ;
z(: ,1) =pinv(C)*y(1:end ,1) ;   %inv((eye(size(A))-A))*(B)*u(1:end ,1) ;

U = zeros(Lu,Ho);

for k=1:(length(yref_opt)/Ly - Ho + 1)

z_Real(1:end,k+1) = A_Real*z_Real(1:end,k) + B_Real*u(1:end,k);
y(1:end,k) = C_Real*z_Real(1:end,k) + D_Real*u(1:end,k);
z(1:end,k+1) = A*z(1:end,k) + B*u(1:end,k);

%y(1:end,k) = C*z(1:end,k) + D*u(1:end,k);
U = QuadSolver(M,P,V,H,f,T,yref_opt(Ly*(k-1)+1:Ly*(k-1)+Ly*Ho),u(1:end,k),z(1:end,k),Ho-1,ymax,ymin,Ly,Lu,umax,umin);
u(1:end,k+1) = U(1:end,2);

k
end
 
%%
close all
figure(1)
subplot(4,1,1)
plot(t,y(1,:),':',t,yref(1,:));
subplot(4,1,2)
plot(t,y(2,:),':',t,yref(2,:));
subplot(4,1,3)
plot(t,y(3,:),':',t,yref(3,:));
subplot(4,1,4)
plot(t,y(4,:),':',t,yref(4,:));

%title("Comparison between yref and the optimized output");
%legend('particle size','solid concentration','sump level' ...
%    ,'circulating load', 'yref particle size','yref solid concentration','yref sump level' ...
%    ,'yref circulating load')

figure(2)
plot(t,u,"x:");
legend('feed rate','mill feed water flow rate','sump dilution water flow rate' ...
    ,'pump speed')
title("Input Control Signal");


figure(3)
plot(t ,y-yref);
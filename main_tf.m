clear all
close all
clc
%% Based on:https://www.sciencedirect.com/science/article/pii/S0892687507001136
%parameters of machine are in the paper

%TODO: check if every tf is correcty written 
%TODO: apply step inputs with the scaled version and compare with paper

%four inputs: 
%feed rate
%mill feed water flow rate
%sump dilution water flow rate
%pump speed

%four controled variables are selected:
%product particle size
%mill solids concentration
%sump level
%circulating load

%%
p11=tf([-0.58],[83 1],'InputDelay',41);
p12=tf([0.97],[125 1],'InputDelay',40)*tf([-1.08*0.97],[195 1],'InputDelay',(56+260));
p13=tf([0.67],[20 1],'InputDelay',8)*tf([-0.67*1.07], [92 1],'InputDelay',214);
p14=tf([0.5], [12 1], 'InputDelay',2);
p21=tf([0.62],[123 1]);
p22=tf([-1.75],[118 1]);
p23=tf([0.51],[81 1], 'InputDelay',87)*tf([1], [182 1]);
p24=tf([0.64],[137 1],'InputDelay',9);
p31=tf([2.61],[110 1], 'InputDelay',45);
p32=tf([9.52],[98 1],'InputDelay',93)*tf([1], [137 1]);
p33=tf([2.83],[128 1],'InputDelay',8)*tf([1], [128 1]);
p34=tf([2.81],[108 1],'InputDelay',8);
p41=tf([1e-3],[150 1],'InputDelay',30)*tf([1],[1 0]);
p42=tf([0.011],[100 1],'InputDelay',30)*tf([1],[1 0]);
p43=tf([0.032],[1 0]);
p44=tf([-0.031],[1 0]);

Pspm=[p11 p12 p13 p14;
      p21 p22 p23 p24;
      p31 p32 p33 p34;
      p41 p42 p43 p44];
Ts = 1;
T = (0:Ts:500)';
figure
bode(Pspm);

figure
impulse(Pspm)
figure
step_resp = step(Pspm, T)


sys_ss = ss(Pspm);
A = sys_ss.A;
B = sys_ss.B;
C = sys_ss.C;
D = sys_ss.D;

%%
close all
d = (dr-0.5)*0.1+1;
step_u = ones(length(T),1) ;
zeros_u = zeros(length(T),1);

y = iddata([step_u,step_u,step_u,step_u],...
           [step_resp(:,1,1,:),step_resp(:,1,2,:),step_resp(:,1,3,:),step_resp(:,1,4,:)],Ts);
y.InputName  = {'step_u';'zeros_u';"zeros_u_1";"zeros_u_2"};
y.OutputName = {'Var5';'Var6';'Var7'; 'Var8'};
dr = rand(1,length(T));
d = ((dr-0.5)*0.1+1)';
step_u = ones(length(T),1) ;
zeros_u = zeros(length(T),1);
y_m = iddata([step_u,d.*step_u,step_u,step_u.*d],...
           [step_resp(:,1,1,:).*d,step_resp(:,1,2,:).*d,step_resp(:,1,3,:).*d,step_resp(:,1,4,:).*d],Ts);
y_m.InputName  = {'step_u';'zeros_u';"zeros_u_1";"zeros_u_2"};
y_m.OutputName = {'Var5';'Var6';'Var7'; 'Var8'};


nx = 3;
sysc = ssest(y,nx , ... 
    'InputName',["step_u" "zeros_u" "zeros_u_1" "zeros_u_2"],...
    'OutputName',["Var5" "Var6" "Var7" "Var8"]);

compare(y_m,y,sysc)

%%
close all
p44=tf([-0.031],[1 0]);
[sys_si, sys_m] = si_task(p44,3,3);
compare(sys_m,sys_si);





clc 
clear all
close all

p11=tf([-0.58],[83 1],'InputDelay',0);
p12=tf([0.97],[125 1],'InputDelay',0)*tf([-1.08*0.97],[195 1],'InputDelay',(0));
p13=tf([0.67],[20 1],'InputDelay',0)*tf([-0.67*1.07], [92 1],'InputDelay',0);
p14=tf([0.5], [12 1], 'InputDelay',0);
p21=tf([0.62],[123 1]);
p22=tf([-1.75],[118 1]);
p23=tf([0.51],[81 1], 'InputDelay',0)*tf([1], [182 1]);
p24=tf([0.64],[137 1],'InputDelay',0);
p31=tf([2.61],[110 1], 'InputDelay',0);
p32=tf([9.52],[98 1],'InputDelay',0)*tf([1], [137 1]);
p33=tf([2.83],[128 1],'InputDelay',0)*tf([1], [128 1]);
p34=tf([2.81],[108 1],'InputDelay',0);
p41=tf([1e-3],[150 1],'InputDelay',0)*tf([1],[1 0]);
p42=tf([0.011],[100 1],'InputDelay',0)*tf([1],[1 0]);
p43=tf([0.032],[1 0]);
p44=tf([-0.031],[1 0]);

sys_real_c=[p11 p12 p13 p14;
      p21 p22 p23 p24;
      p31 p32 p33 p34;
      p41 p42 p43 p44];

p=tf([1], [1  1 1],0.1);

si = siclass();
%si.si_run(p)
sys_id_d(:,1) =p ;
ind = 1;
for i=1:4
    for j=1:4
    [sys_id_d(i,j)]= si.si_run(sys_real_c(i,j));
    i,j
    end
end


sys_id_c = d2c(sys_id_d);
sys_real_d = c2d(sys_real_c,0.1);


bode(sys_id_c, sys_real_c)














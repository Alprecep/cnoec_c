function [sys_id, sys_meas] = si_task(sys,in,out)
%SI_TASK Summary of this function goes here
%   Detailed explanation goes here

Ts = 1;
T = (0:Ts:500)';
dr = rand(1,length(T));
dr1 = rand(1,length(T));
d = ((dr-0.5)*0.1+1)';
d1 = ((dr1-0.5)*0.05+1)';

step_u = ones(length(T),1) ;

sys_resp = step(sys, T);

sys_meas = iddata([step_u.*d1], [sys_resp.*d]);
sys_meas.InputName = {['in' + string(in)]};
sys_meas.OutputName = {['out' + string(out)]};

nx =1:2;
sys_id = ssest(sys_meas,nx,...
    'InputName',['in' + string(in)],...
    'OutputName',['out' + string(out)]);

end


function [delay] = si_delay(sys)
t = 1:1:1000;
response = step(sys, t);
idx = find(response~=0, 1, 'first');
delay = idx-1;

end
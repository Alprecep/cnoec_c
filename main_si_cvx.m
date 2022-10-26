%% si_cvx
clear all
clc
p32=tf([9.52],[98 1],'InputDelay',0)*tf([1], [137 1]); %
p41=tf([1e-3],[150 1],'InputDelay',0)*tf([1],[1 0]); %incorrect

p=tf([1 1], [1 1 1]);

error_list = []

for i=1:3
    for j=1:3
        [a b c e ] = si_cvx(p,i,j);
        error_list(end+1,:) = [e i j];
        
    end
end

varNames = ["error", "na", "nb"];
array2table(error_list, 'VariableNames',varNames)



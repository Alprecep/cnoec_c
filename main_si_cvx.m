clear all
clc
p32=tf([9.52],[98 1],'InputDelay',0)*tf([1], [137 1]); %
p41=tf([1e-3],[150 1],'InputDelay',0)*tf([1],[1 0]); %incorrect

p=tf([9.52],[98 1],'InputDelay',0)*tf([1], [137 1]);

p=tf([1], [1  1 1]);

error_list = []
Ts = 1;
for j=1:5
    for i=j:5
        [a b c d e ] = si_cvx(p,i,j, Ts);
        error_list(end+1,:) = [d, e, i, j, i+j];
        
    end
end

varNames = ["error_t", "error_v", "na", "nb", "compl" ];
error_table = array2table(error_list, 'VariableNames',varNames)
error_table = sortrows(error_table,'compl');

index_min = (find(error_table.error_v == min(error_table.error_v)));
index = find( error_table.error_v < error_table.error_v(index_min)*1.1, 1 );

na = table2array(error_table(index,3));
nb = table2array(error_table(index,4));

na,nb
[system_id, system_real, theta, error_t, error_v] =si_cvx(p,na,nb, Ts);


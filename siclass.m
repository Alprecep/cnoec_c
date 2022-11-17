classdef siclass < handle
    properties
        Tsim = 10;
        Ts = 0.01;
        steps = 1001;
        error_table;
        na = 0;
        nb = 0;
        sys_real;
        sys_id;
        d_in = 40;
        d_out = 40;
        d_val_out = 60;
        d_val_in = 60;

    end
    methods

        function r = u_gen(obj,type)
            if type =="train"
                
                r = ones(1,obj.steps);
                p1 = round(obj.steps/4);
                p2 = 2*p1;
                p3 = 3*p1;
                r(1,1:p1) =  ones(1,p1) + 0.2*square(1:p1, 50);
                r(1,p1+1:p2) =  ones(1,p1) + 0.2*(sin((1:1:p1)/p1*500));
                r(1,p2+1:p3) =  ones(1,p1) + 0.2*rand(1,p1);
            else if type=="validate"
                    r = ones(1,obj.steps) + 0.1*sin(logspace(0 , 1, obj.steps)*5);
                    %r(1,obj.steps) =  ones(1,obj.steps) + 0.2*rand(1,obj.steps);
            else
                r =0;
            end
            end
        end

            function  [system_id] =si_run(obj,sys)
                obj.sys_real = c2d(sys,obj.Ts);

                error_list = [];
                for j=1:4
                    for i=j:4
                        [sys_real, system_id, error_t, error_v] = obj.simopt(sys,i,j);
                        error_list(end+1,:) = [error_t, error_v, i, j, i+j];

                    end
                end

                varNames = ["error_t", "error_v", "na", "nb", "compl"];
                obj.error_table = array2table(error_list, 'VariableNames',varNames);
                obj.error_table = sortrows(obj.error_table,'compl');


                name = 'FPE';
                d = obj.error_table.compl;
                N = obj.Tsim / obj.Ts;

                obj.error_table.(name) = obj.error_table.error_v.*((1+d/N)./(1-d/N));

                index_min = (find(obj.error_table.FPE == min(obj.error_table.FPE)));
                %index = find( obj.error_table.error_v < obj.error_table.error_v(index_min)*1.5, 1 );

                obj.na = table2array(obj.error_table(index_min,3));
                obj.nb = table2array(obj.error_table(index_min,4));

                obj.na;obj.nb;
                [sys_real, system_id, error_t, error_v]=obj.simopt(sys,obj.na,obj.nb);
                obj.sys_id = system_id;

            end

            function [sys, sys_id, error_t,error_v] = simopt(obj,sys,na, nb)
                % timesteps
                T = (0:obj.Ts:obj.Tsim)';
                % disturbances creation
                d = rand(1,length(T))-0.5;
                zvals = obj.u_gen("train");  % input signals
                zvals_d = awgn(zvals,obj.d_in,'measured');

                % generating data for training
                yreal = lsim(sys,(zvals_d), T);
                ymeas = awgn(yreal,obj.d_out,'measured');

                % generating initial guess for the regressors
                x0 = zeros(na+nb,1);
                %% BFGS
                % Initialize solver options
                myoptions               =   myoptimset;
                myoptions.Hessmethod  	=	'BFGS';
                myoptions.gradmethod  	=	'CD';
                myoptions.graddx        =	2^-17;
                myoptions.tolgrad    	=	1e-9;
                myoptions.ls_tkmax      =	1;
                myoptions.ls_beta       =	0.5;
                myoptions.ls_c          =	1e-4;
                myoptions.ls_nitermax   =	100;
                [xopt] = fminunc(@(x) obj.costfcn(na, nb,ymeas, zvals, x, T), x0,myoptions);

                [J, ysim] = obj.costfcn(na, nb, ymeas, zvals, xopt, T);

                % this part is important be aware how to put the parameters of
                % the model. they need to match with regressor phi. and "1" is
                % put since it is a discrete time model -->
                % just write simple time series difference equation convert to
                % z domain tf and you will see the "1"

                theta = xopt;

                sys_id =tf(([theta(na+1:end)])', [1,-1*(theta(1:nb))'], obj.Ts)

                %sys_id= d2c(sys_id); %near zero poles are problem to convert
                %sys;

                %input for validation set
                u_validation =obj.u_gen("validate");
                u_val_d = awgn(u_validation,obj.d_val_in,'measured');

                %validation data generated for real system
                ymeas_v = lsim(sys,(u_val_d), T)  ;
                ymeas_v = ymeas_v(max(na,nb)+1:end);
                ymeas_v = awgn(ymeas_v,obj.d_val_out,'measured');

                %validation for identified data
                yid_v = lsim(sys_id,(u_val_d), T)  ;
                yid_v = yid_v(max(na,nb)+1:end) ;
                yid_v = awgn(yid_v,obj.d_val_out,'measured');
                
                % check the results
                error_t = J;

                res_v = [ymeas_v yid_v];
                error_v = sum(abs(res_v(:,1)-res_v(:,2)));

                sys = c2d(sys,obj.Ts)

            end

            function [J, ysim] = costfcn(obj, na, nb, ymeas, zvals, x, T)
                %%%% SIMULATION %%%%%%%%
                ysim = zeros(size(ymeas));
                J = 0;
                ysim(1:max(na,nb)) = ymeas(1:max(na,nb));

                for ii = max(na,nb)+1:1:size(T,1)
                    ysim(ii) = [ysim(ii-1:-1:ii-na,1)' zvals(ii-1:-1:ii-nb)] * x;
                    % compute the cost
                    J = J + abs((ymeas(ii)-ysim(ii)));
                end
                J = J/size(T,1);
            end

        end
    end

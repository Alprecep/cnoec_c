classdef siclass < handle
    properties
        Tsim = 100;
        Ts = 0.1;
        error_table;
        na = 0;
        nb = 0;
        sys_real;
        sys_id;

    end
    methods

        function [sys_id, sys_real, theta, error_t, error_v] = si_opt(obj,sys,na, nb)
            %SI_TASK Summary of this function goes here
            %   Detailed explanation goes here
            T = (0:obj.Ts:obj.Tsim)';

            % disturbances creation
            dbar = 0;
            dr = rand(1,length(T));
            d = (dr*dbar)';
            zvals = obj.u_gen(T);  % input signals

            d_train_u = d(randperm(length(d)));
            d_train_o = zeros(length(T),1);  % d(randperm(length(d)));
            d_val_u =   d(randperm(length(d)));
            d_val_o =  zeros(length(T),1);  % d(randperm(length(d)));

            % generating data for training
            ymeas = lsim(sys,(zvals+d_train_u'), T) + d_train_o;

            % order of the model
            na_model = na;
            nb_model = nb;

            % assigning the regressor matrix
            PHI=zeros(length(T)-na_model-nb_model,na_model+nb_model);
            for ind=max(na_model,nb_model)+1:length(T)
                PHI(ind,:)      =   [(ymeas(ind-1:-1:ind-na_model))', (zvals(ind-1:-1:ind-nb_model))];
            end
            PHI = PHI(max(na_model,nb_model)+1:end,:);
            ymeas = ymeas(max(na_model,nb_model)+1:end);

            cvx_begin;
            variable t_LASSO(na_model + nb_model,1)
            %minimize norm(t_LASSO,1)
            minimize norm(ymeas-PHI*t_LASSO,inf)
            subject to
            %norm(ymeas-PHI*t_LASSO,inf) <= 5
            cvx_end;

            % optimized variables are removed if they are too small
            theta = t_LASSO;
            for ind=1:length(theta)
                if (norm(theta(ind)) < 1e-10)
                    theta(ind) = 0;
                end
            end

            % this part is important be aware how to put the parameters of
            % the model. they need to match with regressor phi. and "1" is
            % put since it is a discrete time model -->
            % just write simple time series difference equation convert to
            % z domain tf and you will see the "1"

            sys_id =tf(([theta(na_model+1:end)])', [1,-1*(theta(1:na_model))'], obj.Ts);
            %sys_id= d2c(sys_id)
            sys_real = c2d(sys,obj.Ts);

            %input for validation set
            sweep = 2*pi*T/100;
            sweep(1:obj.Tsim*0.3/obj.Ts) = 0;
            u_validation = ones(1,length(T)) + sin(sweep)'*0.1;

            %validation data generated for real system
            ymeas_v = lsim(sys,(u_validation + d_val_u'), T)  ;
            ymeas_v = ymeas_v(max(na_model,nb_model)+1:end);

            %validation for identified data
            yid_v = lsim(sys_id,(u_validation + d_val_u'), T)  ;
            yid_v = yid_v(max(na_model,nb_model)+1:end) ;

            % check the results
            res = [ymeas PHI*theta];
            error_t = sum(abs(res(:,1)-res(:,2)));

            res_v = [ymeas_v yid_v];
            error_v = sum(abs(res_v(:,1)-res_v(:,2)));

        end

        function r = u_gen(obj, T)
            s = length(T);
            r = ones(1,s);
            p1 = round(s/3);
            p2 = 2*p1;
            r(1,1:p1) =  ones(1,p1) + 0.2*square(1:p1, 50);
            r(1,p1+1:p2) =  ones(1,p1) + 0.2*(sin((1:1:p1)/p1*500));
            %r = r*10;
            %r = ones(1,length(T)) + 0.1*sin(logspace(0 , 1, length(T))*5);
             
            %r = ones(1,s);
            % sys identification problem
            %zvals = ones(1,length(T)); % input
            %halfsine = sin(pi*[1:1:10]/10);
            %signal = zeros(1,length(T));
            %signal(1:10) = halfsine;
            %zvals = 50+signal*10;
            %zvals = ones(1,length(T)) + 0.1*sin(logspace(0 , 1, length(T))*5);

            %sweep =log(T');
            %zvals = ones(1,length(T)) + sin(sweep)*0.1; % inp
            %sweep = 2*pi*T/100;
            %sweep(1:Tsim*0.3/Ts) = 0; zvals = ones(1,length(T)) + sin(sweep)'*0.1;
            %zvals = log(T+1)'/10;
            %ymeas =step(p23, T)+d;
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

            varNames = ["error_t", "error_v", "na", "nb", "compl" ];
            obj.error_table = array2table(error_list, 'VariableNames',varNames);
            obj.error_table = sortrows(obj.error_table,'compl');

            index_min = (find(obj.error_table.error_v == min(obj.error_table.error_v)));
            %index = find( obj.error_table.error_v < obj.error_table.error_v(index_min)*1.5, 1 );

            obj.na = table2array(obj.error_table(index_min,3));
            obj.nb = table2array(obj.error_table(index_min,4));

            obj.na,obj.nb;
            [sys_real, system_id, error_t, error_v]=obj.simopt(sys,obj.na,obj.nb);
            obj.sys_id = system_id;

        end

        function [sys, sys_id, error_t,error_v] = simopt(obj,sys,na, nb)
            % timesteps
            T = (0:obj.Ts:obj.Tsim)';
            % disturbances creation
            dbar = 0;
            dr = rand(1,length(T));
            d = (dr*dbar)';
            zvals = obj.u_gen(T);  % input signals

            d_train_u = d(randperm(length(d)));
            d_train_o = zeros(length(T),1);  % d(randperm(length(d)));
            d_val_u =   d(randperm(length(d)));
            d_val_o =  zeros(length(T),1);  % d(randperm(length(d)));

            % generating data for training
            ymeas = lsim(sys,(zvals+d_train_u'), T) + d_train_o;

            % order of the model
            na_model = na;
            nb_model = nb;

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
            
            sys_id =tf(([theta(na_model+1:end)])', [1,-1*(theta(1:na_model))'], obj.Ts)
            sys_real = c2d(sys,obj.Ts)
            
            sys_id= d2c(sys_id);
            sys;

            %input for validation set
            u_validation =ones(1,length(T)) + 0.1*sin(logspace(0 , 1, length(T))*5);
           

            %validation data generated for real system
            ymeas_v = lsim(sys,(u_validation + d_val_u'), T)  ;
            ymeas_v = ymeas_v(max(na_model,nb_model)+1:end);

            %validation for identified data
            yid_v = lsim(sys_id,(u_validation + d_val_u'), T)  ;
            yid_v = yid_v(max(na_model,nb_model)+1:end) ;

            % check the results
            error_t = J;

            res_v = [ymeas_v yid_v];
            error_v = sum(abs(res_v(:,1)-res_v(:,2)));

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
            J = J;
        end

    end
end
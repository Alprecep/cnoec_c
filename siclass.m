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
            dbar = 0.01;
            dr = rand(1,length(T));
            d = (dr*dbar)';
            zvals = obj.u_gen(T);

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
            r(1,1:200) =  ones(1,200) + 0.1*square(1:200, 50);
            r(1,201:500) =  ones(1,300) + 0.1*(sin((1:1:300)/20));
            
            r = ones(1,length(T)) + 0.1*sin(logspace(0 , 1, length(T))*5);
             
            r = ones(1,s);
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
                    [a b c d e ] = obj.si_opt(sys,i,j);
                    error_list(end+1,:) = [d, e, i, j, i+j];

                end
            end

            varNames = ["error_t", "error_v", "na", "nb", "compl" ];
            obj.error_table = array2table(error_list, 'VariableNames',varNames);
            obj.error_table = sortrows(obj.error_table,'compl');

            index_min = (find(obj.error_table.error_v == min(obj.error_table.error_v)));
            index = find( obj.error_table.error_v < obj.error_table.error_v(index_min)*1.2, 1 );

            obj.na = table2array(obj.error_table(index,3));
            obj.nb = table2array(obj.error_table(index,4));

            obj.na,obj.nb;
            [system_id, sys_real, theta, error_t, error_v]=obj.si_opt(sys,obj.na,obj.nb);
            obj.sys_id = system_id;

        end


    end
end
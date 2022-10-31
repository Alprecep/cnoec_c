function [sys_id, sys_real, theta, error_t, error_v] = si_cvx(sys,na, nb, Ts)
    %SI_TASK Summary of this function goes here
    %   Detailed explanation goes here
    %Ts = 0.1;
    Tsim=1000;
    T = (0:Ts:Tsim)';
    
    % disturbances creation
    dbar = 0.001;
    dr = rand(1,length(T));
    d = (dr*dbar)';
    
    % sys identification problem
    %zvals = ones(1,length(T)); % input
    %halfsine = sin(pi*[1:1:10]/10);
    %signal = zeros(1,length(T));
    %signal(1:10) = halfsine; 
    %zvals = 50+signal*10;
    %zvals = ones(1,length(T)) + 0.1*sin(logspace(0 , 1, length(T))*5);
    
    zvals = ones(1,length(T)) + 0.01*square(T/5, 50)';
    zvals(1,1:100) =  ones(1,100);
    zvals(1,100:499) =  ones(1,400) + 0.01*(sin((1:1:400)/20));
    
    %sweep =log(T');
    %zvals = ones(1,length(T)) + sin(sweep)*0.1; % inp
    %sweep = 2*pi*T/100;
    %sweep(1:Tsim*0.3/Ts) = 0; zvals = ones(1,length(T)) + sin(sweep)'*0.1;
    %zvals = log(T+1)'/10;
    %ymeas =step(p23, T)+d;                     % measures    
    ymeas = lsim(sys,(zvals+ + d'), T) ;
    
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
    
    % optimized variables
    theta = t_LASSO;
    for ind=1:length(theta)
        if (norm(theta(ind)) < 1e-10)
            theta(ind) = 0;
        end
    end
    
    sys_id =tf(([theta(na_model+1:end)])', [1,-1*(theta(1:na_model))'], Ts);
    
    
    ymeas_id = lsim(sys_id,zvals, T) + d;
    display("Here is the result")
    sys_id
    %sys_id= d2c(sys_id)
    sys_real = c2d(sys,Ts)
    
    sweep = 2*pi*T/100;
    sweep(1:Tsim*0.3/Ts) = 0; 
    zvals = ones(1,length(T)) + sin(sweep)'*0.1;
    
    ymeas_v = lsim(sys,(zvals+ + d'), T) ;
    ymeas_v = ymeas_v(max(na_model,nb_model)+1:end);
    
    yid_v = lsim(sys_id,(zvals+ + d'), T) ;
    yid_v = yid_v(max(na_model,nb_model)+1:end);
    
    % check the results
    res = [ymeas_v PHI*theta];
    error_t = sum(abs(res(:,1)-res(:,2)));
    
    res_v = [ymeas_v yid_v];
    error_v = sum(abs(res_v(:,1)-res_v(:,2)));
    
    end
    
    
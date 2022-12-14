function U = QuadSolver(M,P,V,H,f,T,yref,U0,Z0,N,ymax,ymin,Ly,Lu,umax,umin)
dy1 = [];
dy2 = [];
O = [];
du1 = [];
du2 = [];

for k = 1:N
    y = yref((k-1)*Ly + 1:k*Ly,1) - yref((k)*Ly + 1:(k+1)*Ly,1);
    O = [O ; y];
end 

% Construction of the constraint vector for the control signal
for i=1:N
    du1 = [du1 ; (umax - U0)];
    du2 = [du2 ; (U0 - umin)];
end 
du = [du1 ; du2];

Z = P*Z0 + V*U0 + O;
for i=1:(N+1)
    dy1 = [dy1 ; (yref((i-1)*Ly+1:i*Ly,1) - ymin)];
    dy2 = [dy2 ; (ymax - yref((i-1)*Ly+1:i*Ly,1))];
end 

d = [dy1 ;  dy2 ; du];

sol = quadprog(H,f,T,d,M,Z);

%sol = quadprog(H,f,[],[],M,Z);

U = zeros(Lu,N);
DeltaU = sol((N+1)*Ly + 1:end);
U(1:end,1) =  DeltaU(1:Lu,1) + U0;
for k=2:N
    U(1:end,k) =  DeltaU((k-1)*Lu+1:k*Lu,1) + U(1:end,k-1);
end

U = [U0 U];

end
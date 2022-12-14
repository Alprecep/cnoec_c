function [M,P,V,H,f,T] = ConstantMatrices(A,B,C,D,Q,R,N,Ly,Lu)

Iy = eye(Ly);
YlinesUcolumns = zeros(Ly,Lu);
L = [];
M = zeros(N*length(Iy), length(Iy)*(N + 1) + N*Lu);
P = [];
V = [];

% M matrix construction
for k = 1:N

    % Construction of the k-th line part containing [0 0 ... 0 I -I 0 ... 0]
    i = 1;
    while(i <= N+1)
        switch i
            case k
                switch k
                    case 1
                        L = [Iy -Iy];
                    case N+1
                        L = [L 0*Iy];
                    otherwise
                        L = [L Iy -Iy];
                end
                i=i+1;
            otherwise
                L = [L 0*Iy];
        end
        i=i+1;
    end
    M(1+(k-1)*Ly:k*Ly,1:(N+1)*Ly) = L;
    L = [];

    % Construction of the k-th line part containing [CA^(k-2)B CA^(k-3)B ... CA^(1)B CB D 0 ... 0]
    switch k
        case 1
            for j=1:N
                if(j == 1)
                    L = D;
                else
                    L = [L YlinesUcolumns];
                end
            end
        otherwise
            for j=1:(k-1)
                L = [L C*A^(k-2 - (j - 1))*B];
            end
            L = [L D];
            for j=(k+1):N
                L = [L YlinesUcolumns];
            end
    end
    M(1+(k-1)*Ly:k*Ly,(N+1)*Ly + 1:end) = -L;
    L = [];
end

% P matrix Construction
for j=1:N
    P = [P ; C*(A^(j) - A^(j-1))];
end

% V matrix Construction
for j=1:N
    V = [V ; C*A^(j-1)*B];
end

% Construction of H and f
H = zeros((N+1)*Ly + N*Lu);
QLinesRColumns = zeros(Ly,Lu);
RLinesQColumns = zeros(Lu,Ly);


% Q sector
for k=1:(N+1)

    for i=1:(N+1)
        switch i
            case k
                L = [L Q];
            otherwise
                L = [L 0*Q];
        end
    end
    for i=1:N
        L = [L QLinesRColumns];
    end

    H(1+(k-1)*Ly:k*Ly,1:end) = L;
    L = [];
end

% R sector
for k=1:(N)

    for i=1:(N+1)
        L = [L RLinesQColumns];
    end

    for i=1:N
        switch i
            case k
                L = [L R];
            otherwise
                L = [L 0*R];
        end
    end

    H(1+(k-1)*Lu+(N+1)*Ly:k*Lu+(N+1)*Ly,1:end) = L;
    L = [];
end

f = zeros((N+1)*Ly + N*Lu,1);

% Ty+ matrice construction
Ty1 = zeros((N+1)*Ly);
for i=1:(N+1)
    for k=1:(N+1)
        switch k
            case i
                L = [L Iy];
            otherwise
                L = [L 0*Iy];
        end
    end
    Ty1((i-1)*Ly+1:i*Ly,1:end) = L;
    L = [];
end

% Ty matrice construction
Ty = [Ty1; -Ty1];

% Tu+ matrice construction
Tu1 = zeros(N*Lu);
Iu = eye(Lu);
for i=1:(N)
    for k=1:(N)
       if(k <= i)
           L = [L Iu];
       else 
           L = [L 0*Iu];
       end 
    end
    Tu1((i-1)*Lu+1:i*Lu,1:end) = L;
    L = [];
end

% Ty matrice construction
Tu = [Tu1; -Tu1];

% T matrice construction
O1 = zeros(2*(N+1)*Ly,N*Lu); 
O2 = zeros(2*N*Lu,(N+1)*Ly);
T = [Ty O1; O2 Tu];
end



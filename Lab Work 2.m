a = 4;
b = 9;
if (a > b)
    error('a should be less than b');
end

var = 15;
n_arr = [4,6,8,10,12];
error_list = zeros(size(n_arr, 2));

for sizeIt = 1:size(n_arr,2)
    size_n = n_arr(sizeIt);
    
    [A, f] = FillMatrixAndFCol(var, size_n, a, b);
    
    n = size(A,1);
    
    S = zeros(size(A));
    D = zeros(size(A));
    
    D(1,1) = sign(A(1,1));
    S(1,1) = sqrt(A(1,1)); 
    
    for i = 2:n
        D(i,i) = sign(A(i,i) - sum((S(1:i-1, i).^2)' * D(1:i-1,1:i-1)) ); 
        S(i,i) = sqrt(abs( A(i,i) - sum((S(1:i-1, i).^2)' * D(1:i-1,1:i-1)) ) ); 
    end
    
    for i = 1:n
        for j = i + 1:n
            TS = 0;
            for k = 1:i-1
                TS = TS + S(k,i)*D(k,k)*S(k,j);
            end
            S(i,j) = ( A(i,j) - TS )/(S(i,i)*D(i,i));
        end
    end
    
    x = FindSolution(n, f, S, D);
    
    error_list(sizeIt) = norm (f - A * x) / norm(f);
    
end

plot(n_arr, error_list);

function [m, f] = FillMatrixAndFCol(var, size_n, a, b)
    m = zeros(size_n);
    for i = 1:size_n 
        for j = 1:size_n
            m(i, j) = 1 / (i + j + var);
        end
    end
    
    f = zeros(size_n, 1);
    for i = 1:size_n 
        f(i) = a + (b - a) * rand; 
    end
end

function [x] = FindSolution(n, f, S, D)
    y = FindY(n, f, S, D);

    x = zeros(n,1);
    x(n) = y(n)/S(n,n);
    for i = n-1:-1:1
        TS = 0;
        for j = i+1:n
            TS = TS + S(i,j)*x(j);
        end
        x(i) = (y(i) - TS)/S(i,i);
    end
end

function [y] = FindY(n, f, S, D)
    y = zeros(n,1);
    y(1) = f(1) / (S(1,1)*D(1,1));
    for i=2:n
        TS = 0;
        for j=1:i-1
            TS = TS + S(j,i)*y(j)*D(j,j);
        end
        y(i) = (f(i)-TS)/(S(i,i)*D(i,i));
    end
end
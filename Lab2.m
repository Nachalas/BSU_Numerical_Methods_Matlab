% Lab 2
% Variant 10
% a(i,j) = 1 / (i + j + v)
% f = a + (b - a) * rand, where a,b are random && a<b 

while 1
    a = input('Enter the variable a: ');
    b = input('Enter the variable b: ');
    if (a < b)
        break
    end
end
variant = 10;

n_options = [4,6,8,10,12];
% idx = randperm(length(n_options),1);
discrepancy_list = zeros(size(n_options, 2)); % Список отн. норм невязок для каждой матрицы, которая будет рассмотрена

for size_iteration = 1:size(n_options,2)
    temp_n = n_options(size_iteration);
    A = zeros(temp_n);
    f = zeros(temp_n, 1);
    for i = 1:temp_n % Initialization of a matrix and f column
        f(i) = a + (b - a) * rand; % Initializing f column
        for j = 1:temp_n
            A(i, j) = 1 / (i + j + variant);
        end
    end
    
    x = SquareRoot(A,f);
    
    discrepancy_list(size_iteration) = norm (f - A * x) / norm(f);
    
end

Plot(discrepancy_list, n_options)

function [x] = SquareRoot(A,f)
  S = zeros(size(A));
  D = zeros(size(A));
  n = size(A,1);
  
  D(1,1) = sign(A(1,1)); % sum equals 0 when i = 0 in (5)
  S(1,1) = sqrt(A(1,1)); % sum equals 0 when i = 0 in (6)
  
  for i = 2:n
    % diagonal elements of D and S
    D(i,i) = sign(A(i,i) - sum((S(1:i-1, i).^2)' * D(1:i-1,1:i-1)) ); %(5)
    S(i,i) = sqrt(abs( A(i,i) - sum((S(1:i-1, i).^2)' * D(1:i-1,1:i-1)) ) ); %(6)
  end
  
  for i = 1:n
    for j = i + 1:n
      temp_sum = 0;
      for k = 1:i-1
        temp_sum = temp_sum + S(k,i)*D(k,k)*S(k,j);
      end
      S(i,j) = ( A(i,j) - temp_sum )/(S(i,i)*D(i,i));
    end    
  end
  
  y = find_y(f,S,D,n);
  x = find_x(y, S,n);
  end

function Plot(discrepancy_list, n)
    plot(n, discrepancy_list);
    xlabel("Matrix size");
    ylabel("Relative norm of discrepancy");
end

function [y] = find_y(f, S, D, n)
  y = zeros(n,1);
  y(1) = f(1) / (S(1,1)*D(1,1)); 
  for i=2:n
    temp_sum = 0;
    for j=1:i-1
      temp_sum = temp_sum + S(j,i)*y(j)*D(j,j);
    end
    y(i) = (f(i)-temp_sum)/(S(i,i)*D(i,i));
  end
end

function [x] = find_x(y,S,n)
  x = zeros(n,1);
  x(n) = y(n)/S(n,n);
  for i = n-1:-1:1
    temp_sum = 0;
    for j = i+1:n
      temp_sum = temp_sum + S(i,j)*x(j);
    end
    x(i) = (y(i) - temp_sum)/S(i,i);
  end
end
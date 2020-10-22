eps = 1e-3;
var = 10;

a = 2;
b = 5;
if (a > b)
    error('a is not less than b');
end

n = 10;

[A, f] = Fill_A_and_F(n, a, b, var);

det_of_A = norm(A);
taus = [1/(2*det_of_A), 1/(4*det_of_A), 1/(8*det_of_A)];

for i = 1:size(taus, 2)
   H = eye(n) - taus(i) * A; 
   phi = taus(i) * f;
   
   q = norm(H);
   xs = phi;
   xn = H * xs + phi;
   
   iteration_epsilon = (q / (1 - q)) * norm(xn - xs);
   eps_list = [iteration_epsilon];
   
   while(iteration_epsilon >= eps)
      xs = xn;
      xn = H * xs + phi;
      iteration_epsilon = (q / (1 - q)) * norm(xn - xs);
      eps_list = AppendEps(eps_list, iteration_epsilon);
   end
   
   Plot(eps_list);
   ConsoleOutput(xn, i);
   
end

legend('tau = 1/(2*||A||)', 'tau = 1/(4*||A||)', 'tau = 1/(6*||A||)');

function [eps_list] = AppendEps(eps_list, iteration_epsilon)
    eps_list = [eps_list; iteration_epsilon];
end

function [] = ConsoleOutput(xn, i)
    disp(strcat('Solution column number:', string(i)));
    disp(xn);
end

function [] = Plot(eps_list)
   hold on;
   plot(eps_list);
end

function [A, f] = Fill_A_and_F(n, a, b, var)
    A = zeros(n);
    f = zeros(n, 1);

    for i = 1:n
        f(i) = a + (b - a) * rand;
        for j = 1:n
            if(i == j)
                A(i,j) = 100 + var;
            else
                A(i, j) = 1 / (i + j + var);
            end
        end
    end
end
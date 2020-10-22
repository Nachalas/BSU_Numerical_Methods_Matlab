a = input('Enter a: ');
b = input('Enter b: ');
if (b < a)
    error('b should be not less than a!');
end

v = 15;
n = 10;
eps = 1e-3;

A = FillMatrix(n, v);
f = FillFColumn(n, a, b);

normOfA = norm(A);
tauVariants = [1/(2*normOfA), 1/(4*normOfA), 1/(8*normOfA)];

for i = 1:size(tauVariants, 2)
   H = eye(n) - tauVariants(i) * A; 
   phi = tauVariants(i) * f;
   
   q = norm(H);
   xs = phi;
   xn = H * xs + phi;
   
   to_plot = CountEpsilons(q, xs, xn, eps, H, phi);
   
   disp(strcat('row number:', string(i)));
   PlotAndOutput(to_plot, xn);
   
end

legend('tau = 1/(2*||A||)', 'tau = 1/(4*||A||)', 'tau = 1/(6*||A||)');

function [to_plot] = CountEpsilons(q, xs, xn, eps, H, phi)
    temp_e = (q / (1 - q)) * norm(xn - xs);
    to_plot = [temp_e];
   
   while(temp_e >= eps)
      xs = xn;
      xn = H * xs + phi;
      temp_e = (q / (1 - q)) * norm(xn - xs);
      to_plot = [to_plot; temp_e];
   end
end

function [] = PlotAndOutput(to_plot, xn)
   hold on;
   plot(to_plot);
   disp(xn);
end

function [A] = FillMatrix(n, v)
    A = zeros(n);
    for i = 1:n
    for j = 1:n
        if(i == j)
            A(i,j) = 100 + v;
        else
            A(i, j) = 1 / (i + j + v);
        end
    end
    end
end

function [f] = FillFColumn(n, a, b)
    f = zeros(n, 1);
    for i = 1:n
           f(i) = a + (b - a) * rand; 
    end
end
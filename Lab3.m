% Lab 2
% Variant 10
% a(i,j) = 1 / (i + j + v) (i != j), a(i,j) = 100 + v (i == j)
% f = rand in range (a,b), where a,b are random && a<b 

while 1
    a = input('Enter the variable a: ');
    b = input('Enter the variable b: ');
    if (a < b)
        break
    end
end

variant = 10;
n = 10;
eps = 0.001;

A = zeros(n);
f = zeros(n, 1);
for i = 1:n
    f(i) = a + (b - a) * rand;
    for j = 1:n
        if(i == j)
            A(i,j) = 100 + variant;
        else
            A(i, j) = 1 / (i + j + variant);
        end
    end
end

A_norm = norm(A);
tau_list = [1/(2*A_norm), 1/(4*A_norm), 1/(8*A_norm)];
Identity_matr = eye(n);

for i = 1:size(tau_list, 2)
   H = Identity_matr - tau_list(i) * A; 
   phi = tau_list(i) * f;
   
   q = norm(H);
   xs = phi;
   xn = H * xs + phi;
   
   curr_eps = (q / (1 - q)) * norm(xn - xs);
   to_plot = [curr_eps];
   
   while(curr_eps >= eps)
      xs = xn;
      xn = H * xs + phi;
      curr_eps = (q / (1 - q)) * norm(xn - xs);
      to_plot = [to_plot; curr_eps];
   end
   
   hold on;
   plot(to_plot);
   disp(xn);
   
end

legend('tau = 1/(2*||A||)', 'tau = 1/(4*||A||)', 'tau = 1/(6*||A||)');

%disp(A)
%disp(f)
%disp(Identity_matr);
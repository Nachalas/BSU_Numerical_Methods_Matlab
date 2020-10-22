clear;

n = 7; N = n^2;
eps = 1.e-5;
A = gallery('poisson', n);

f = rand(N,1);

max_iteration_count = 1000;

B = diag(diag(A));
[Pogr_1, k_1] = FindEpsJacobyOrZeidel(B, A, f, N, max_iteration_count, eps);

B = triu(A);
[Pogr_2, k_2] = FindEpsJacobyOrZeidel(B, A, f, N, max_iteration_count, eps);

[Pogr_3, k_3] = FindEpsTriag(A, f, N, max_iteration_count, eps);

Plot(k_1, Pogr_1, k_2, Pogr_2, k_3, Pogr_3);

function [data, k] = FindEpsJacobyOrZeidel(B, A, f, N, max_iteration_count, eps)
    InvB = inv(B);
    D = InvB * A;
    g = InvB * f;
    x = zeros(N, 1);
    r = D * x - g;
    eps_t = norm(r) / norm(g);    
    k = 0;
    while (eps_t > eps && k < max_iteration_count)
           A_mult_r = A * r;
           tau = (r' * r)/((A_mult_r') * r);
           x = x - tau*r;
           r = D * x - g;
           k = k + 1;
           eps_t = norm(r)/norm(g);
           data(k) = eps_t;
    end
end

function [data, k] = FindEpsTriag(A, f, N, max_iteration_count, eps) 
R1 = GetRMatrix(A, N, 1);
R2 = GetRMatrix(A, N, 0);
w = 0.5;
B = (eye(N) + w * R1) * (eye(N) + w * R2);
    InvB = inv(B);
    D = InvB * A;
    g = InvB * f;
    x = zeros(N, 1);
    r = D * x - g;
    eps_t = norm(r) / norm(g);    
    k = 0;
    while (eps_t > eps && k < max_iteration_count)
           A_mult_r = A * r;
           tau = (r' * r)/((A_mult_r') * r);
           x = x - tau*r;
           r = D * x - g;
           k = k + 1;
           eps_t = norm(r)/norm(g);
           data(k) = eps_t;
    end
end

function [R] = GetRMatrix(A, N, bool)
R = zeros(N);
for k = 1:N
   for m = 1:N
       if (k > m && bool == 1)
          R(k, m) = A(k, m);
       elseif(k < m && bool == 0)
          R(k, m) = A(k, m);
       elseif (k == m)
          R(k, m) = A(k, m) / 2;
       end
   end
end
end

function [] = Plot(k_1, Pogr_1, k_2, Pogr_2, k_3, Pogr_3)
semilogy(1:k_1, Pogr_1, 1:k_2, Pogr_2, 1:k_3, Pogr_3)
legend('Jacoby','Zeidel','Triag');
grid
end
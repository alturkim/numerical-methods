A = [5 4 1 1;4 5 1 1;1 1 4 2;1 1 2 4];
x = ones(length(A),1);
% using Power Method
[k1, x1, r1] = power_method(A,x, 2);
disp("result of my implementation of Power Method");
disp("number of iterations:");
disp(k1);
disp("Resulted Eigenvalue");
disp(r1);
disp("-----------------------------------------------------------------");
% using accelerated Power Method
[k1, x1, r1] = accelerated_power_method(A,x, 2);
disp("result of my implementation of Accelerated Power Method");
disp("number of iterations:");
disp(k1);
disp("Resulted Eigenvalue");
disp(r1);
disp("-----------------------------------------------------------------");
% using Inverse Power Method
[k2, x2, r2] = inverse_power_method(A,x,2);
disp("result of my implementation of Inverse Power Method");
disp("number of iterations:");
disp(k2);
disp("Resulted Eigenvalue");
disp(r2);
disp("-----------------------------------------------------------------");
% using Accelerated Inverse Power Method
[k2, x2, r2] = accelerated_inverse_power_method(A,x,2);
disp("result of my implementation of Accelerated Inverse Power Method");
disp("number of iterations:");
disp(k2);
disp("Resulted Eigenvalue");
disp(r2);
disp("-----------------------------------------------------------------");
% using Shifted Power Method
[k3, x3, r3] = shifted_power_method(A,x,3, 2);
disp("result of my implementation of Shifted Power Method, mu=3");
disp("number of iterations:");
disp(k3);
disp("Resulted Eigenvalue");
disp(r3);
disp("-----------------------------------------------------------------");
% using Accelerated Shifted Power Method
[k3, x3, r3] = accelerated_shifted_power_method(A,x,3, 2);
disp("result of my implementation of Accelerated Shifted Power Method, mu=3");
disp("number of iterations:");
disp(k3);
disp("Resulted Eigenvalue");
disp(r3);
disp("-----------------------------------------------------------------");
% using Shifted Inverse Power Method
[k4, x4, r4] = shifted_inverse_power_method(A,x, (r1+r2)/2, 2);
disp("result of my implementation of Shifted Inverse Power Method");
disp("number of iterations:");
disp(k4);
disp("Resulted Eigenvalue");
disp(r4);
disp("-----------------------------------------------------------------");
% using Accelerated Shifted Inverse Power Method
[k4, x4, r4] = accelerated_shifted_inverse_power_method(A,x, (r1+r2)/2, 2);
disp("result of my implementation of Accelerated Shifted Inverse Power Method");
disp("number of iterations:");
disp(k4);
disp("Resulted Eigenvalue");
disp(r4);
disp("-----------------------------------------------------------------");
% verifying the result using matlab's function eig()
[V,D] = eig(A);
disp("Results from matlab's function eig()");
disp(V);
disp(D);
disp("__________________________________________________________________");



function [k, x, r] = power_method(A, x, linear_func_idx)
epsilon = 1e-15;
kmax = 10000;
r = 0;
for k=1:kmax
   r_prev = r;
   y = A*x;
   r = y(linear_func_idx)/x(linear_func_idx);
   x = y/norm(y,inf);
   if norm(r-r_prev) < epsilon && r ~=0
      return;
   end
end
end

function [k, x, s] = accelerated_power_method(A, x, linear_func_idx)
epsilon = 1e-15;
kmax = 10000;
s = 0;
for k=1:kmax
   if k>=3
      r_prev_prev = r_prev;
   end
   if k>=2
      r_prev = r;
   end
   y = A*x;
   r = y(linear_func_idx)/x(linear_func_idx);
   x = y/norm(y,inf);
   if k>=3
      s_prev = s;
      s = r - (r-r_prev)^2/(r-(2*r_prev)+r_prev_prev);
      if norm(s-s_prev) < epsilon  && s ~=0
          return;
       end
   end
end
end



function [k, x, r] = inverse_power_method(A, x, linear_func_idx)
    epsilon = 1e-15;
    kmax = 10000;
    r = 0;
    for k=1:kmax
       r_prev = r;
       y = linsolve(A,x);
       r = x(linear_func_idx)/y(linear_func_idx);
       x = y/norm(y,inf);
       if norm(r-r_prev) < epsilon && r ~=0
          return;
       end
    end
end

function [k, x, s] = accelerated_inverse_power_method(A, x, linear_func_idx)
epsilon = 1e-15;
kmax = 10000;
s = 0;
for k=1:kmax
   if k>=3
      r_prev_prev = r_prev;
   end
   if k>=2
      r_prev = r;
   end
   y = linsolve(A,x);
   r = x(linear_func_idx)/y(linear_func_idx);
   x = y/norm(y,inf);
   if k>=3
      s_prev = s;
      s = r - (r-r_prev)^2/(r-(2*r_prev)+r_prev_prev);
      if norm(s-s_prev) < epsilon  && s ~=0
          return;
       end
   end
end
end


function [k, x, r] = shifted_power_method(A, x, mu, linear_func_idx)
A = A - mu*eye(length(A));
[k, x, r] = power_method(A, x, linear_func_idx);
r = r+mu;
end


function [k, x, r] = shifted_inverse_power_method(A, x, mu, linear_func_idx)
    A = A - mu*eye(length(A));
    [k, x, r] = inverse_power_method(A, x, linear_func_idx);
    r = r + mu;
end

function [k, x, r] = accelerated_shifted_power_method(A, x, mu, linear_func_idx)
A = A - mu*eye(length(A));
[k, x, r] = accelerated_power_method(A, x, linear_func_idx);
r = r+mu;
end


function [k, x, r] = accelerated_shifted_inverse_power_method(A, x, mu, linear_func_idx)
A = A - mu*eye(length(A));
[k, x, r] = accelerated_inverse_power_method(A, x, linear_func_idx);
r = r + mu;
end



A = [7 1 -1 2; 1 8 0 -2; -1 0 4 -1; 2 -2 -1 6];
b = [3; -5; 4; -3];
x = [0; 0; 0; 0];
[k, x] = jacobi(A, b, x);
disp("Using Jacobi");
disp(k);
disp(x);

function [k, x] = jacobi(A, b, x)
    kmax = 2000;
    epsilon = 5e-5;
    delta = 1e-15;
    n = length(A);
    for k=1:kmax
       y = x;
       for i=1:n
          sum = b(i);
          diag = A(i,i);
          if abs(diag) < delta
            disp("diagonal element too small");
            return;
          end
          for j=1:n
             if j ~= i
                sum = sum - (A(i,j)*y(j)); 
             end
          end
          x(i) = sum/diag;
       end
       if norm(x-y,inf) < epsilon
           return
       end
    end
    disp("maximum iterations reached");
end



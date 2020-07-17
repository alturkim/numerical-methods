A = [7 1 -1 2; 1 8 0 -2; -1 0 4 -1; 2 -2 -1 6];
b = [3; -5; 4; -3];
x = [0 0 0 0];
[k, x] = gauss_seidel(A, b, x);
disp("Using Gauss-Seidel");
disp(k);
disp(x);


function [k, x] = gauss_seidel(A, b, x)
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
          for j=1:i-1
            sum = sum - (A(i,j)*x(j)); 
          end
          for j=i+1:n
            sum = sum - (A(i,j)*x(j)); 
          end
          x(i) = sum/diag;
       end
       if norm(x-y,inf) < epsilon
           return
       end
    end
    disp("maximum iterations reached");
end

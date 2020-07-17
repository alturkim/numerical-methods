A = [4 6 10; 6 25 19; 10 19 62]; 
disp("using matrix A");
% calling my function:
disp("output of my function");
Cholesky(A)
disp("output of matlab function");
% calling the function in matlab library to check
chol(A,'lower')

B = [4 6 10; 6 13 19;10 19 62];
disp("using matrix B");
% calling my function:
disp("output of my function");
Cholesky(B)
disp("output of matlab function");
% calling the function in matlab library to check
chol(B,'lower')


function [L] = Cholesky(A)
    n = length(A);
    L = zeros(n,n);
    for k=1:n
        sum = 0;
        for s=1:k-1
            sum = sum + L(k,s)^2;
        end
        L(k,k) = sqrt(A(k,k) - sum);

        sum = 0;
        for i=k+1:n
           for s=1:k-1
              sum = sum + (L(i,s) * L(k,s)); 
           end
           L(i,k) = (A(i,k) - sum)/L(k,k); 
        end
    end
end

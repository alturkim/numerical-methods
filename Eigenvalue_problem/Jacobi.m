A = [1 4; 4 1];
disp("Test Case 1: A=");
disp(A);
disp("Output of My implementation of Jacobi's Method");
[V, D] = jacobi(A)
disp("Output of Matlab's eig()");
[V, D] = eig(A)
disp("------------------------------------------------------------------");

A = [2 0 3;0 4 1;3 1 1];
disp("Test Case 2: A=");
disp(A);
disp("Output of My implementation of Jacobi's Method");
[V, D] = jacobi(A)
disp("Output of Matlab's eig()");
[V, D] = eig(A)


function [V, D] = jacobi(A)
    D = A;
    n = length(A);
    V = eye(n);
    kmax = 10000000;
    epsilon = 1e-5;

    for k=1:kmax        
        max_off_diag = 0;
        p = 0;
        q = 0;
        for i=1:n
            for j=1:n
                if i ~= j
                    if abs(D(i,j)) >= max_off_diag
                        max_off_diag = abs(D(i,j));
                        p = i;
                        q = j;
                    end
                end
            end
        end
        
        theta = (1/2)* acotd((D(q,q)-D(p,p))/(2*D(p,q)));
        S = eye(n);
        S(p,p) = cos(theta);
        S(q,q) = cos(theta);
        S(p,q) = sin(theta);
        S(q,p) = -sin(theta);
 
        D = S' * D * S;
        V = V*S;
        
        % calculate off diagonal sum of squares to check for convergence
        off_diag_sum = 0;
        for i=1:n
            for j=1:n
                if i~=j
                    off_diag_sum = off_diag_sum + D(i,j)^2;
                end
            end
        end
        
        % check for convergence
        if sqrt(off_diag_sum) < epsilon
            return;
        end
        
    end

end
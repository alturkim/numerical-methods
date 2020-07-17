A = [5 7 6 5;7 10 8 7;6 8 10 9;5 7 9 10];
n = length(A);

% Case 1: Doolittle
L = zeros(n,n);
U = zeros(n,n);
for i=1:n
   L(i,i) = 1; 
end
[L, U] = factorize(A, L, U);
disp("case 1: output");
disp("L");
disp(L);
disp("U");
disp(U);
disp("verify LU = A");
x = L*U;
disp(x);

% --------------------------------------------------------
% Case 2: Crout
L = zeros(n,n);
U = zeros(n,n);
for i=1:n
   U(i,i) = 1; 
end
[L, U] = factorize(A, L, U);
disp("case 2: output");
disp("L");
disp(L);
disp("U");
disp(U);
disp("verify LU = A");
x = L*U;
disp(x);
% --------------------------------------------------------



function [L, U] = factorize(A, L, U)
    n = length(A);
    
    if L(1,1)~=0
        U(1,1) = A(1,1)/L(1,1);
    elseif U(1,1)~=0
        L(1,1) = A(1,1)/U(1,1);
    end
  
    for i=1:n
       L(i,1) = A(i,1)/U(1,1);
       U(1,i) = A(1,i)/L(1,1);
    end
    
    
    for k=2:n
        sum = 0;
        for m=1:k-1
            sum = sum + (L(k,m)*U(m,k));
        end
        if L(k,k)~=0
            U(k,k) = (A(k,k)-sum)/L(k,k);
        elseif U(1,1)~=0
            L(k,k) = (A(k,k)-sum)/U(k,k);
        end
        
        for i=k:n
           sum = 0;
           for m=1:k-1
             sum = sum + L(i,m)*U(m,k);
           end
           L(i,k) = (1/U(k,k)) * (A(i,k) - sum);
        end
    

    for j=k:n
       sum = 0;
       for m=1:k-1
         sum = sum + L(k,m)*U(m,j);
       end
       U(k,j) = (1/L(k,k)) * (A(k,j) - sum);
    end
    end
end

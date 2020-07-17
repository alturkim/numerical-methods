
A = [0.0001 -5.0300 5.8090 7.8320;
     2.2660 1.9950 1.2120 8.0080;
     8.8500 5.6810 4.5520 1.3020;
     6.7750 -2.2530 2.9080 3.9700];
 b = [9.5740; 7.2190; 5.7300; 6.2910];
 solve(A, b, "naive")
 solve(A, b, "spp")

function [x] = solve(A, b, mode)
    switch mode
        case "naive"
            x = naive_gauss(A,b);
        case "pp"
            x = pp_gauss(A,b);
        case "spp"
           x = spp_gauss(A,b);
    end
end

function [x] = naive_gauss(A, b)
    n = length(A);
    m = length(b(1,:));
    x = zeros(n,m);
    for k=1:n-1
       for i=k+1:n
           xmult = A(i,k)/A(k,k);
           A(i,k) = xmult;
           for j=k+1:n
               A(i,j) = A(i,j) - (xmult*A(k,j));
           end
           b(i,:) = b(i,:) - (xmult*b(k,:));
       end
    end
    x(n,:) = b(n,:)/A(n,n);
    for i=n-1:-1:1
        sum = b(i,:);
        for j=i+1:n
            sum = sum - (A(i,j)*x(j,:));
        end
        x(i,:) = sum/A(i,i);
    end
end

function [x] = pp_gauss(A, b)
    n = length(A);
    m = length(b(1,:));
    x = zeros(n,m);
    L = 1:n;
    
    for k=1:n-1
       rmax = 0;
       for i=k:n
          r = abs(A(L(i),k));
          if r>rmax
            rmax = r;
            j = i;
          end
       end
       temp = L(j);
       L(j) = L(k);
       L(k) = temp;
       for i=k+1:n
           xmult = A(L(i),k)/A(L(k),k);
           A(L(i),k) = xmult;
           for j=k+1:n
               A(L(i),j) = A(L(i),j) - (xmult*A(L(k),j));
           end
           b(L(i),:) = b(L(i),:) - (xmult*b(L(k),:));
       end
    end
    x(n,:) = b(L(n),:)/A(L(n),n);
    for i=n-1:-1:1
        sum = b(L(i),:);
        for j=i+1:n
            sum = sum - (A(L(i),j)*x(j,:));
        end
        x(i,:) = sum/A(L(i),i);
    end
end

function [x] = spp_gauss(A, b)
    n = length(A);
    m = length(b(1,:));
    x = zeros(n,m);
    L = 1:n;
    s = double.empty(0,n);
    
    for i=1:n
        smax = 0;
        for j=1:n
            smax = max(smax, abs(A(i,j)));
        end
        s(i) = smax;
    end
    
    for k=1:n-1
       rmax = 0;
       for i=k:n
          r = abs(A(L(i),k)/s(L(i)));
          if r>rmax
            rmax = r;
            j = i;
          end
       end
       temp = L(j);
       L(j) = L(k);
       L(k) = temp;
       for i=k+1:n
           xmult = A(L(i),k)/A(L(k),k);
           A(L(i),k) = xmult;
           for j=k+1:n
               A(L(i),j) = A(L(i),j) - (xmult*A(L(k),j));
           end
           b(L(i),:) = b(L(i),:) - (xmult*b(L(k),:));
       end
    end
    x(n,:) = b(L(n),:)/A(L(n),n);
    for i=n-1:-1:1
        sum = b(L(i),:);
        for j=i+1:n
            sum = sum - (A(L(i),j)*x(j,:));
        end
        x(i,:) = sum/A(L(i),i);
    end
end

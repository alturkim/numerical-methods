t = 0;
x = 0;
b = 1;
h = 2^-7;
n = (b-t)/h;
[X, T] = RK4(@f, t, x, h, n);

function val = f(t, x)
    if t == 0 && x == 0
        val = 0;
    else
        val = (x/t) + t*sec(x/t);
    end
end

function [X, T] = RK4(f, t, x, h, n)
    exact_sol = @(t) t*asin(t); 
    X = double.empty(n+1,0);
    T = double.empty(n+1,0);
    a = t;
    X(1) = x;
    T(1) = t;
    for j=1:n
        K1 = h*f(t, x);
        K2 = h*f(t+(0.5)*h, x+(0.5)*K1);
        K3 = h*f(t+(0.5)*h, x+(0.5)*K2);
        K4 = h*f(t+h, x+K3);
        
        x = x+(1/6)*(K1+ 2*K2+ 2*K3+ K4);
        t = a + j*h;
        X(j+1) = x;
        T(j+1) = t;
        exact = exact_sol(t);
        fprintf("t= %f\tx= %f\texact= %f\n",t,x,exact);
    end
end


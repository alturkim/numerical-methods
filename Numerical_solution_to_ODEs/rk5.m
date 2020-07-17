t = 0;
x = 0;
b = 1;
h = 2^-7;
n = (b-t)/h;
[X, T] = RK5(@f, t, x, h, n);

function val = f(t, x)
    if t == 0 && x == 0
        val = 0;
    else
        val = (x/t) + t*sec(x/t);
    end
end


function [X, T] = RK5(f, t, x, h, m)
    exact_sol = @(t) t*asin(t); 
    X = double.empty(m+1,0);
    T = double.empty(m+1,0);
    a = t;
    X(1) = x;
    T(1) = t;
    for j=1:m
        K1 = h*f(t, x);
        K2 = h*f(t+(0.5)*h, x+(0.5)*K1);
        K3 = h*f(t+(0.5)*h, x+(0.25)*K1+(0.25)*K2);
        K4 = h*f(t+h, x-K2+2*K3);
        K5 = h*f(t+(2/3)*h, x+(7/27)*K1+(10/27)*K2+(1/27)*K4);
        K6 = h*f(t+(1/5)*h, x+(28/625)*K1-(1/5)*K2+(546/625)*K3+(54/625)*K4-(378/625)*K5);
        
        x = x + (1/24)*K1 + (5/48)*K4 + (27/56)*K5 + (125/336)*K6;
        t = a + j*h;
        X(j+1) = x;
        T(j+1) = t;
        exact = exact_sol(t);
        fprintf("t= %f\tx= %f\texact= %f\n",t,x,exact);
    end
end

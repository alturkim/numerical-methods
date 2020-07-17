f = @(x) x+x^2;
exact_sol = @(t) exp(t)/(16-exp(t));
x = exp(1)/(16-exp(1));
a = 1;
b = 2.77;
h = 1/100;
n = (b-a)/h;
[T, X] = Taylor(a, x, n, h);
function [T, X] = Taylor(a, x, n, h)
    exact_sol = @(t) exp(t)/(16-exp(t));
    X = double.empty(n+1,0);
    T = double.empty(n+1,0);
    t = a;
    T(1) = t;
    X(1) = x;
    for k=1:n   
       x_p1 = x + x^2;
       x_p2 = x_p1 + 2*x*x_p1;
       x_p3 = x_p2 + 2*x_p1^2 + 2*x*x_p2;
       x_p4 = x_p3 + 6*x_p1*x_p2 + 2*x*x_p3;
       x_p5 = x_p4 + 6*x_p2^2 + 8*x_p1*x_p3 + 2*x*x_p4;
       
       x = x + h*(x_p1 + (h/2)*(x_p2 + (h/3)*(x_p3+(h/4)*(x_p4+(h/5)*(x_p5)))));
       t = a + k*h;
       exact = exact_sol(t);
       fprintf("t= %f\tx= %10.10f\texact= %10.10f\n",t,x,exact);
       T(k+1) = t;
       X(k+1) = x;
    end
        
    
end

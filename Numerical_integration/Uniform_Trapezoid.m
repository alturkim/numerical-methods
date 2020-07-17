f = @(x) exp((-x^2));
a = 0;
b = 1;
n = 60;
sum = Trapezoid_uniform(f, a, b, n);
fprintf("integral of exp((-x^2)) in [0,1] with 60 uniform intervals = %f\n",sum);

function sum = Trapezoid_uniform(f, a, b, n)
    h = (b-a)/n;
    sum = 0.5*(f(a)+f(b));
    for i=1:n-1
       x = a+i*h;
       sum = sum + f(x);
    end
    sum = sum*h;
end

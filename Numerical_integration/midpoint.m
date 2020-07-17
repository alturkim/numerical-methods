f = @(x) exp(-x)*cos(x);
a = 0;
b = 2*pi;
n = 32;

disp("midpoint rule result:");
sum = midpoint_rule(f, a, b, n);
disp(sum);


function sum = midpoint_rule(f, a, b, n)
    h = (b-a)/n;
    sum = 0;
    for i=0:n-1
       x = a+i*h;
       sum = sum + f(x+0.5*h);
    end
    sum = sum*h;
end








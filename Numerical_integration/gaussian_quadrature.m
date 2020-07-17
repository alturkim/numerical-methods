f = @(x) (x^-1)*sin(x);
a = 0;
b = 1;
sum = three_point_gauss(f,a,b,1000);
fprintf("integrating (x^-1)*sin(x) from 0 to 1 with 1000 subinterval\nresult=%10.10f\n",sum);
function sum = three_point_gauss(f,a,b, n)
    h = (b-a)/n;
    sum = 0;
    for i=0:n-1
        c1 = a+i*h;
        c2 = c1+h;
        g = @(t) 0.5*(c2-c1)*t + 0.5*(c2+c1);
        sum = sum + 0.5*(c2-c1)*((5/9)*f(g(-sqrt(3/5))) + (8/9)*f(g(0)) +(5/9)*f(g(sqrt(3/5))));
    end
end








f = @(x) cos(2*x)/exp(x);
a = 0;
b = 2*pi;
n = 120;
true_sol = (1/5)*(1-exp(-2*pi));

disp("Simpson's 1/3  result:");
sum = simpson_third(f,a,b, n);
error = abs(true_sol -sum);
fprintf("Result=%f\t |error|=%e\n", sum, error);
disp("------------------------------------------------------------------");


function sum = simpson_third(f,a,b, n)
    h = (b-a)/n;
    temp = 0;
    for i=1:0.5*n-1
        x = a+(2*i-1)*h;
        temp = temp + 2*f(x);
        x = a + 2*i*h;
        temp = temp + f(x);
    end
    sum = temp*2*h/3;
    sum = sum + (h/3)*(f(a)+f(b)) + (4*h/3)*f(b-h);
end









format long
disp("Using equally spaced nodes");
x = linspace(-1,1,14);
y = atan(x);
t = linspace(-1,1,100);
disp("coefficients in the Newton form of the polynomial");
a = Coef(x,y);
for i=1:length(a)
    fprintf("a_%d=%.10f\n",i-1,a(i));
end

maxerror = 0;
for i=1:length(t)
    val = Eval(x,a,t(i));
    y = atan(t(i));
    if abs(y-val)>= maxerror
        maxerror = abs(y-val);
    end
end
fprintf("Maximum absolute error |f(x) - p(x)| over 100 points= %e\n",maxerror)

disp("-----------------------------------------------------------------------");
disp("Using Chebyshev nodes");
x = double.empty(14,0);
for i=0:13
    x(i+1) = cos((2*i+1)*pi/30);
end
y = atan(x);
t = linspace(-1,1,100);
disp("coefficients in the Newton form of the polynomial");
a = Coef(x,y);
for i=1:length(a)
    fprintf("a_%d=%.10f\n",i-1,a(i));
end

maxerror = 0;
for i=1:length(t)
    val = Eval(x,a,t(i));
    y = atan(t(i));
    if abs(y-val)>= maxerror
        maxerror = abs(y-val);
    end
end
fprintf("Maximum absolute error |f(x) - p(x)| over 100 points= %e\n",maxerror)


function a = Coef(x, y)
    n = length(x);
    a = zeros(n,1);
    for i=1:n
        a(i) = y(i);
    end
    for j=2:n
        for i=n:-1:j
            a(i) = (a(i) - a(i-1))/(x(i) - x(i-j+1));
        end
    end       
end

function val = Eval(x,a, t)
    n = length(x);
    val = a(n);
    for i=n-1:-1:1
        val = val*(t-x(i)) + a(i);
    end
end

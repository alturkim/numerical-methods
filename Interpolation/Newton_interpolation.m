format long
x = linspace(1,6,11);
y = atan(x);
t = linspace(0,8,33);
disp("coefficients in the Newton form of the polynomial");
a = Coef(x,y);
for i=1:length(a)
    fprintf("a_%d=%.10f\n",i-1,a(i));
end

for i=1:length(t)
    val = Eval(x,a,t(i));
    y = atan(t(i));
    fprintf("arctan(%f)=%f\tp(%f)=%f\tarctan(x)-p(x)=%e\n",t(i),y,t(i),val, y-val);
end

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

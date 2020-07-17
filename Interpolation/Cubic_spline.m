t = [-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0   1  2  3  4  5 6 7 8 9 10 11];
y = [ -6 -5 -6 -5 -4 -5 -3 -4 -3 -1 -3 -2 -1 -2 -3 -1 1 0 1 2 1   2];

z = Spline3_Coef(t,y);

Xs = linspace(-10,11,100);
Ys = double.empty(length(Xs),0);
for i=1:length(Xs)
    Ys(i) = Spline3_Eval(t, y, z, Xs(i));
end
plot(Xs, Ys);
grid on
grid minor

    




function z = Spline3_Coef(t,y)
    n = length(t);
    z = double.empty(n,0);
    u = double.empty(n-1,0);
    v = double.empty(n-1,0);
    h = double.empty(n-1,0);
    b = double.empty(n-1,0);
    for i=1:n-1
        h(i) = t(i+1) - t(i);
        b(i) = (y(i+1) - y(i))/h(i);
    end
    u(1) = 2*(h(1)+h(2));
    v(1) = 6*(b(2)-b(1));
    for i=2:n-1
        u(i) = 2*(h(i) + h(i-1)) - (h(i-1)^2/u(i-1));
        v(i) = 6*(b(i) - b(i-1)) - (h(i-1)*v(i-1)/u(i-1));
    end
    z(n) = 0;
    for i=n-1:-1:2
        z(i) = (v(i)-h(i)*z(i+1))/u(i);
    end
    z(1) = 0;
end

function val = Spline3_Eval(t, y, z, x)
    n = length(t);
    for i=n-1:-1:1
        if x-t(i) >= 0
            break;
        end
    end
    h = t(i+1)-t(i);
    temp = (z(i)/2) + (x-t(i))*(z(i+1)-z(i))/(6*h);
    temp = -(h/6)*(z(i+1)+2*z(i)) + (y(i+1) - y(i))/h + temp*(x-t(i));
    val = y(i) + (x-t(i))*temp;
end

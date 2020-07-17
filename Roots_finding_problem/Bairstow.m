disp("Solving x^4 - 5x^3 + 10x^2 - 10x + 4 = 0");
a = [4, -10, 10, -5, 1];
[xs] = BairstowMethod(a)
disp("-------------------------------------------------------------");
disp("Solving x^4 + 2x^3 - 11x^2 + 8x -60 = 0");
a = [-60, 8, -11, 2, 1];
[xs] = BairstowMethod(a)
disp("-------------------------------------------------------------");
disp("Solving 6x^5 + 11x^4 - 33x^3 - 33x^2 + 11x + 6");
a = [6, 11, -33, -33, 11, 6];
[xs] = BairstowMethod(a)
disp("-------------------------------------------------------------");

function [xs] = BairstowMethod(a)
    n = length(a);
    xs = zeros(n-1,1);
    u = - a(n-1)/a(n);
    v = - a(n-2)/a(n);
    
    for step=1:2:4
        [u, v, a] = update(u,v,a);
        xs(step) = (u+sqrt(u^2+4*v))/2;
        xs(step+1) = (u-sqrt(u^2+4*v))/2;
        a = a(3:end);
        if length(a) == 3
           xs(end-1) = (-a(2)+sqrt(a(2)^2-4*a(3)*a(1)))/(2*a(3));
           xs(end) = (-a(2)-sqrt(a(2)^2-4*a(3)*a(1)))/(2*a(3));
           return;
        elseif length(a) == 2
            xs(end) = -a(1)/a(2);
            return;
        end
    end
    
    
end

function [u, v, b] = update(u,v, a)
    n = length(a);
    for step=1:100
        bs = zeros(n,1);
        cs = zeros(n,1);
        b(n) = a(n);
        b(n-1) = a(n-1)+b(n)*u;
        c(n) = b(n);
        c(n-1) = b(n-1) + c(n)*u;
        for i=n-2:-1:1
           b(i) = a(i) + b(i+1)*u + b(i+2)*v;
           c(i) = b(i) + c(i+1)*u + c(i+2)*v;
        end
        oldu = u;
        oldv = v;
        u = ((b(1)*c(4)-b(2)*c(3))/(c(3)^2-c(2)*c(4))) + u;
        v = ((b(2)*c(2)-b(1)*c(3))/(c(3)^2-c(2)*c(4))) + v;
        if abs(oldu - u) < 1e-15
            return;
        end
    end
end
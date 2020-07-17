% f = @(x) x^3 + 2*x^2 +10*x - 20;
f = @(x) -20+(x*(10+(x*(2+x))));
% f_prime = @(x) 3*x^2 + 4*x +10;
f_prime = @(x) 10+(x*(4+(3*x)));
x = 2;
[n, x, fx] = Newton(f, f_prime, x, 10, 0.5e-5, 0.5e-5)
function [n, x, fx] = Newton(f, f_prime, x, nmax, epsilon, delta)
    fx = f(x);
    for n=1:nmax
        n
        fp = f_prime(x)
        if abs(fp) < delta
            disp("small derivative");
            return;
        end
        d = fx/fp;
        x = x - d
        fx = f(x)
        if abs(d) < epsilon
            disp("Newton's method converged");
            return;
        end
    end
end

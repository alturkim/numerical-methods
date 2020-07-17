f = @(x) x^3 - x^2 - 2*x + 1;
a = -0.5;
b = 0.5   ;
nmax = 100;
epsilon = 0.5e-6    ;
[n, c, fc, error] = Bisection(f, a, b, nmax, epsilon)

function [n, c, fc, error] =  Bisection(f, a, b, nmax, epsilon)
    fa = f(a);
    fb = f(b);
    if sign(fa) == sign(fb)
        disp("function has same signs at a and b");
        return;
    end
    error = b - a;
    for n=0:nmax
        error = error/2;
        c = a + error;
        fc = f(c);
        
        if abs(error) < epsilon
            disp("The method converged");
            return;
        end
        
        if sign(fa) ~= sign(fc)
            b = c;
            fb = fc;
        else
            a = c; 
            fa = fc;
        end
    end

end

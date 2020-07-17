f = @(x) exp(x) - 3*x^2;
x1 = -1;
x2 = 1;
nmax = 100;
[n, a, fa, ierr] = Secant(f, x1, x2, 1e-54, 1e-54, nmax);
fprintf("n=%i\nx=%f\nfx=%f\n",n,a,fa);

function [n, x1, fa, ierr] = Secant(f, x1, x2, epsilon, delta, nmax)
    ierr = 0;
    fa = f(x1);
    fb = f(x2);
    if abs(fa) > abs(fb)
        temp = x1;
        x1 = x2;
        x2 = temp;
        temp = fa;
        fa = fb;
        fb = temp;
    end
    for n=2:nmax
        if abs(fa) > abs(fb)
            temp = x1;
            x1 = x2;
            x2 = temp;
            temp = fa;
            fa = fb;
            fb = temp;
        end
        d = (x2-x1)/(fb-fa);
        x2 = x1;
        fb = fa;
        d = d*fa;
        
        if abs(d) < epsilon
            disp("Method Stopped:The method converged");
            ierr = 1;
            return;
        end
        x1 = x1 - d;
        fa = f(x1);
        if abs(fa) < delta
            ierr = 1;
            disp("Method Stopped:Function value is too small");
            return;
        end
    end
       
        
end

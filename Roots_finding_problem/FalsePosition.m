% finding the root of f(x)
disp("Finding a root for f(x)");
f = @(x) x^3 + 3*x -1;
a = 0;
b = 1;
nmax = 50;
epsilon = 1e-54;

disp("Using False Position");
[n, c, fc, error] = FalsePosition(f, a, b, nmax, epsilon)
disp("___________________________________________________________________");



function [n, c, fc, error] =  FalsePosition(f, a, b, nmax, epsilon)
    fa = f(a);
    fb = f(b);
    if sign(fa) == sign(fb)
        disp("function has same signs at a and b");
        return;
    end
    error = b - a;
    for n=0:nmax
        error = error/2;
        c = (a*f(b) - b*f(a))/(f(b) - f(a));
        fc = f(c);
        
        if abs(error) < epsilon
            disp("The method converges");
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



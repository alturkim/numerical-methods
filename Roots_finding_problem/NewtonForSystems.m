% part a
disp("part a");
X = [1; 1];
F = @(X)[X(1)^2+X(2)^2-25;
         X(1)^2-X(2)-2];
J = @(X)[2*X(1),2*X(2);
         2*X(1), -1];
[n, X, FX] = NewtonForSystems(F, J, X, 100)


function [n, X, FX] = NewtonForSystems(F, J, X, nmax)
    FX = F(X);
    for n=1:nmax
        JX = J(X);
        H = linsolve(JX,FX);
        X = X - H; 
        FX = F(X);
    end
end

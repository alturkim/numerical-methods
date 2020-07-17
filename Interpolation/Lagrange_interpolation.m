x = linspace(0,1.6875, 20);
y = sin(x);
t = linspace(0,1.6875, 40);
disp("coefficients in the Lagrange form of the polynomial");
a = Coef(x);
for i=1:length(a)
    fprintf("a_%d=%.10f\n",i-1,a(i));
end

for i=1:length(t)
    val = Eval(x,y,a,t(i));
    sinx = sin(t(i));
    fprintf("sin(%f)\t-\tp(%f)\t=\t%e\n",t(i),t(i), sinx-val);
end



function a = Coef(x)
    n = length(x);
    a = ones(1,n);
    for i=1:n
       for j=1:n
           if i~=j
              a(i) = a(i) * 1/(x(i)-x(j));
           end
       end
    end
end

function val = Eval(x,y, a, t)
    n = length(x);
    val = 0;
    for i=1:n
        temp = a(i);
        for j=1:n
            if i ~= j
                temp = temp * (t-x(j));
            end
        end
        val = val + y(i) * temp;
    end
end

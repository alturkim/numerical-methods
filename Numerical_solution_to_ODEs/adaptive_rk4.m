t1 = 0;
x1 = 0;
tb = 10;
h = 0.01;
itmax = 1000;
ep_max = 1e-5;
ep_min = 1e-8;
h_min = 1e-6;
h_max = 1;
f = @(x,t) 3 + 5*sin(t) + 0.2*x;
fprintf("solving x'(t)= 3 + 5*sin(t) + 0.2*x with x(0)=0\n")
% calling RK4 with one iteration to get a second initial value
[x2, t2] = RK4(f, t1, x1, h);

[X, T] = adaptive(f, t1, x1, t2, x2, h, tb, itmax, ep_max, ep_min, h_min, h_max);

for i=1:length(X)
    fprintf("t=%10.15f\t x=%f\n",T(i),X(i));
end

function [X, T] = adaptive(f, t1, x1, t2, x2, h,tb, itmax, ep_max, ep_min, h_min, h_max)
    delta = 0.5e-5;
    iflag = 1;
    k=2;
    itmax = itmax + 2;
    X = double.empty(itmax, 0);
    T = double.empty(itmax, 0);
    X(1) = x1;
    X(2) = x2;
    T(1) = t1;
    T(2) = t2;
    
    while k<=itmax
        k=k+1;
        if abs(h) < h_min
            h = sign(h)*h_min;
        end
        if abs(h) > h_max
            h = sign(h)*h_max;
        end
        d = abs(tb-T(k-1));
        
        if d<=abs(h)
            iflag = 0;
            if d <= delta*max(abs(tb),abs(T(k-1)))
                break;
            end
            h = sign(h)*d;
        end
        x_d = X(k-1) + (h/2)*(3*f(T(k-1),X(k-1))-f(T(k-2),X(k-2)));
        x = X(k-1) + (h/2)*(f(T(k-1)+h, x_d)+f(T(k-1),X(k-1)));
        t = T(k-1) + h;
        ep = (1/6)*abs(x-x_d);
        
        if iflag ==0
            X(k) = x;
            T(k) = t;
            break
        end
        if ep < ep_min
            h = 2*h;
            X(k) = x;
            T(k) = t;
        elseif ep > ep_max
            h = h/2;
            k = k-1;
        else
            X(k) = x;
            T(k) = t;
        end
    end
end

function [x, t] = RK4(f, t, x, h)
    a = t;
    K1 = h*f(t, x);
    K2 = h*f(t+(0.5)*h, x+(0.5)*K1);
    K3 = h*f(t+(0.5)*h, x+(0.5)*K2);
    K4 = h*f(t+h, x+K3);
    x = x+(1/6)*(K1+ 2*K2+ 2*K3+ K4);
    t = a + h;
end

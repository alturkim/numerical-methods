format long

f = @(x) (1+x)^-1;
a=0;
b=1;
n=10;
R = Romberg(f,a,b,n);
I = log(1+b)-log(1+a);
diff = I-R(n,n);
fprintf("The error in calculating integral of (1+x)^-1 in [%d,%d] with n=%d is:%e\n",a,b,n,diff);
disp("---------------------------------------------------------------------");

f = @(x) exp(x);
a=0;
b=1;
n=10;
R = Romberg(f,a,b,n);
I = exp(b)-exp(a);
diff = I-R(n,n);
fprintf("The error in calculating integral of e^x in [%d,%d] with n=%d is:%e\n",a,b,n,diff);
disp("---------------------------------------------------------------------");

f = @(x) (1+x^2)^-1;
a=0;
b=1;
n=10;
R = Romberg(f,a,b,n);
I = atan(b) - atan(a);
diff = I-R(n,n);
fprintf("The error in calculating integral of (1+x^2)^-1 in [%d,%d] with n=%d is:%e\n",a,b,n,diff);
disp("---------------------------------------------------------------------");

function R = Romberg(f, a, b, n)
    R = zeros(n,n);
    h = b-a;
    R(1,1) = (h/2)*(f(a) + f(b));
    
    for i=2:n
        h = h/2;
        sum = 0;
        for k=1:2:(2^(i-1))-1
            sum = sum + f(a+k*h);
        end
        R(i,1) = (1/2)*R(i-1,1) + (sum*h);
        for j=2:i
            R(i,j) = R(i,j-1) + (R(i,j-1) - R(i-1,j-1))/(4^(j-1) - 1);
        end
    end
end

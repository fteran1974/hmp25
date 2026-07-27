function[A,B]=fiedler(coeff,v)
% constructs the coefficients, A,B of the Fiedler linearization \la A+ B
% associated with a polynomial whose coefficients are given by coeff and
% whose PCIS is v
% coeff is a matrix of size nxn(p+1), where p is the degree of the
% polynomial

[n r]= size(coeff);
p = r/n-1;

A = eye(n*p);
A(1:n,1:n) = coeff(:,n*p+1:end);
if p == 0
    B = eye(n);
    p = 1;
else
    if v(1)==1
        B=[-coeff(:,n+1:2*n),eye(n);-coeff(:,1:n),zeros(n)];
    else
        B=[-coeff(:,n+1:2*n),-coeff(:,1:n);eye(n),zeros(n)];
    end
    for i=1:p-2
        if v(i+1)==1
            B=[-coeff(:,n*(i+1)+1:n*(i+2)),eye(n),zeros(n,n*i);
                B(:,1:n),zeros(n*(i+1),n),B(:,n+1:end)];
        else
            B=[-coeff(:,n*(i+1)+1:n*(i+2)),B(1:n,:);
                eye(n),zeros(n,n*(i+1));zeros(n*i,n),B(n+1:end,:)];
        end
    end
end
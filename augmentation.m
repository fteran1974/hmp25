function[evals,invcond,rvec,lvec,alpha,beta] = augmentation(coeff,k,delta)
% PEP solved by augmentation in [Hochstenbach-Mehl-Plestenjak,2025]
% coeff=[A0 ... Ap] are the coefficients of the matrix poly P
% n-k is the rank of P
% delta is the threshold (default eps^(1/2)
% evals is the array of eigenvalues
% invcond is the reciprocal of the condition number
% rvec,lvec are the corresponding left and right eigenvectors
% v is a vector of length p-1 which corresponds to the
% consecution-inversion sequence of the Fiedler linearization which is used
% to solve the associated GEP
[n,r]=size(coeff);
p = r/n-1; % degree of P
U=orth(rand(n,k));V=orth(rand(n,k));
Q1=rand(k,k*(p+1));Q2=rand(k,k*(p+1));
pa=zeros((n+k),(p+1)*(n+k));
for i=1:p+1
    pa(:,(n+k)*(i-1)+1:(n+k)*i)=[coeff(:,n*(i-1)+1:n*i),U*Q1(:,k*(i-1)+1:k*i);Q2(:,k*(i-1)+1:k*i)*V',zeros(k)];
end
v=int32(randi([0, 1], [1, p-1]));
[X,Y,E,~,sp,~,~,~,~]=scalfiedlereigconbis(pa,v);
%[X,Y,E,~,sp,~,~,~,~]=scalfiedlereigcon(coeff,v);
alpha=zeros([1,(n+k)*p]);
for i=1:(n+k)*p
    alpha(i)=norm(X(n+1:end,i));
end
beta=zeros([1,(n+k)*p]);
for i=1:(n+k)*p
    beta(i)=norm(Y(n+1:end,i));
end
[~,ind]=find(max(alpha,beta)<delta);

%nu = zeros(p+1,1);
%for i = 1:p+1, nu(i) = norm(coeff(:,n*(i-1)+1:n*i),'fro'); end
%s = zeros(n*p,1);  
% Compute condition numbers as in [Hochstenbach-Mehl-Plestenjak,2025].
for j = 1:(n+k)*p
    la = (E(j).^(0:p-1));
    Der = coeff(:,n+1:2*n);
    for k = 2:p
        Der = Der+k*la(k)*coeff(:,n*k+1:n*(k+1));
    end
    nabp = norm((E(j).^(0:p)));%norm((E(j).^(0:p)).*nu');
    % condition number for the polynomial:
    sp(j) = nabp*norm(Y(1:n,j))*norm(X(1:n,j))/abs(Y(1:n,j)'*Der*X(1:n,j));
    %nab = norm([a(j)*norm(A),b(j)*norm(B)]);
end

evals=E(ind);
invcond=1./sp(ind);
rvec=X(:,ind);
lvec=Y(:,ind);
alpha=alpha(ind);
beta=beta(ind);

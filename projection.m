function[evals,infinite,gamma,gammainf,alpha,beta] = projection(coeff,k,delta,delta2,delta3,delta4)
% PEP solved by projection in [Hochstenbach-Mehl-Plestenjak,2025]
% coeff=[A0 ... Ap] are the coefficients of the matrix poly P
% n-k is the rank of P
% delta is the threshold (default eps^(1/2)
% evals is the array of eigenvalues
% invcond is the reciprocal of the condition number
% gamma is the reciprocal of the condition number computed as in
% [HMP25]
% rvec,lvec are the corresponding left and right eigenvectors
% v is a vector of length p-1 which corresponds to the
% consecution-inversion sequence of the Fiedler linearization which is used
% to solve the associated GEP
[n,r]=size(coeff);
p = r/n-1;
U=orth(randn(n));
V=orth(randn(n));
U1=U(:,1:n-k);V1=V(:,1:n-k);
U2=U(:,n-k+1:end);V2=V(:,n-k+1:end);
% Construct the projection matrix poly
% Compute the eigenvalues and normalized right and left eigenvectors of the
% projected matrix poly
coeffp=zeros(n-k,(n-k)*(p+1));
for i=1:p+1
coeffp(:,(i-1)*(n-k)+1:i*(n-k))=U1'*coeff(:,(i-1)*n+1:i*n)*V1;
end
v=int32(randi([0, 1], [1, p-1]));
[X,Y,E]=scalfiedlereigconbis(coeffp,v);
% compute alphas and betas
for i=1:p*(n-k)
    pol=coeff(:,1:n);
   for j=1:p
        pol=pol+coeff(:,n*j+1:n*(j+1))*E(i)^j;
    end
    alpha(i)=norm(U2.'*pol*V1*X(:,i));
    beta(i)=norm(Y(:,i).'*U1'*pol*V2);%V2.'*pol'*U1*Y(:,i));
end
[~,ind]=find(max(alpha,beta)<delta);
% Compute condition numbers as in [Hochstenbach-Mehl-Plestenjak,2025].
for j = 1:p*(n-k)
    la = (E(j).^(0:p-1));
    Der = coeffp(:,(n-k)+1:2*(n-k));
    for s = 2:p
        Der = Der+s*la(s)*coeffp(:,(n-k)*s+1:(n-k)*(s+1));
    end
    nabp = norm((E(j).^(0:p)));%norm((E(j).^(0:p)).*nu');
    sp(j) = nabp*norm(Y(:,j))*norm(X(:,j))/abs(Y(:,j).'*Der*X(:,j));
end
gamma=1./sp;
% Compute the gap as in [HMP25]
for i=1:(n-k)*2
    dist=abs(E-E(i));
    dist=dist(dist>0);
    gap(i)=min(dist)/sqrt(1+abs(E(i))^2);
end
ind2=[];
for i=1:(n-k)*2
if gamma(i)<delta2
    ind2=[ind2,i];
elseif gamma(i)<delta3 && gap(i)>delta4
    ind2=[ind2,i];
end
end
in=isnan(gamma);
[~,in]=find(in==1);
ind2=union(ind2,in);
infinite=E(ind2);
gammainf=gamma(ind2);
[~,ind3]=find(max(alpha,beta)>delta);
ind4=union(ind2,ind3);
ind=setdiff(1:2*(n-k),ind4);
evals=E(ind);
invcond=1./sp(ind);
rvec=X(:,ind);
lvec=Y(:,ind);
alpha=alpha(ind);
beta=beta(ind);

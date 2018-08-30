function [nabla]=GetNabla(G,gam,X)

if nargin==2
	X=Xinfd(G,gam);
end

[A,B1,B2,C1,C2,D11,D12,D21,D22]=GetSS(G);
[n,p,q,l,m]=Getsz(G);

B=[B1 B2];
Cb=[C1;zeros(l,n)];
Db=[D11 D12;eye(l) zeros(l,m)];
J=[eye(p) zeros(p,l); zeros(l,p) eye(l)*-gam^2];

Rb=Db'*J*Db+B'*X*B;
	
R1=Rb(1:l,1:l);
R2=Rb(1:l,l+1:l+m)';
R3=Rb(l+1:l+m,l+1:l+m);

nabla=R1-R2'*inv(R3)*R2;

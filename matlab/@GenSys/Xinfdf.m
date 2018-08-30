function [X,nIndefNabla]=Xinfdf(G,gam,Xf,N)
% [X,nIndefNabla]=Xinfdf(G,gam,Xf,N)
%
% Computes solution to the Hinf Discrete Riccati Difference Equation asscoiated
% with P. Ref: Limebeer and Green 1995
%
% N:  horizon length
% Xf: terminal cost

[A,B1,B2,C1,C2,D11,D12,D21,D22]=GetSS(G);
[n,p,q,l,m]=Getsz(G);

B=[B1 B2];
Cb=[C1;zeros(l,n)];
Db=[D11 D12;eye(l) zeros(l,m)];
J=[eye(p) zeros(p,l); zeros(l,p) eye(l)*-gam^2];

X=Xf;
nIndefNabla=0;
for i=N:-1:1
	Rb=Db'*J*Db+B'*X*B;
	Lb=Cb'*J*Db+A'*X*B;
	X=Cb'*J*Cb+A'*X*A-(Lb*inv(Rb))*Lb';

	R1=Rb(1:l,1:l);
	R2=Rb(1:l,l+1:l+m)';
	R3=Rb(l+1:l+m,l+1:l+m);

	nabla=R1-R2'*inv(R3)*R2;

	if max(eig(nabla))>0
		nIndefNabla=nIndefNabla+1;
	end
end






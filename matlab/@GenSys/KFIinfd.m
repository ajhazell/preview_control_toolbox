function [K,F2,F0]=KFIinfd(P,X,gam)
% K=KFIinfd(P,X)
%
% Generate full information, discrete-time, Hinf-suboptimal  controller gain

[A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(P);
[n,p,q,l,m]=Getsz(P);

B=[B1 B2];
Cb=[C1;zeros(l,n)];
Db=[D11 D12;eye(l) zeros(l,m)];
J=[eye(p) zeros(p,l); zeros(l,p) eye(l)*-gam^2];

R=Db'*J*Db+B'*X*B;
L=Db'*J*Cb+B'*X*A;


R1=R(1:l,1:l);
R2=R(1:l,l+1:l+m)';
R3=R(l+1:l+m,l+1:l+m);
L2=L(l+1:l+m,:);
invR3=inv(R3);
F2=-invR3*L2;
F0=-invR3*R2;
K=[F2 F0];
function [A,B,R,L,Q,Cb,Db,J]=GetSSb(sys,gam)
% Get parameters for full information Hinfinity Riccati equation
[A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys);
[n,p,q,l,m]=Getsz(sys);

B=[B1 B2];
Cb=[C1;zeros(l,n)];
Db=[D11 D12;eye(l) zeros(l,m)];
J=[eye(p) zeros(p,l); zeros(l,p) eye(l)*-gam^2];
R=Db'*J*Db;
L=Cb'*J*Db;
Q=Cb'*J*Cb;



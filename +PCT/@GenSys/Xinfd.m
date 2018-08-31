function [X,F,R,L]=Xinfd(P,gam)
% [X,F,R,L,info]=Xinfd(P,gam)
% Compute solution to the Hinfinity Full Information Riccati equation
% associated with P
%
% X:    DARE solution
% F:    feedback gain 
% L:    closed-loop eigenvalues
% R:    Rbar in Thesis

[A,B1,B2,C1,C2,D11,D12,D21,D22]=GetSS(P);
[n,p,q,l,m]=Getsz(P);

%E'XE=A'XA-(A'XB+S)(B'XB+R)^-1(A'XB+S)'+Q
%
%A=A
%E=I
%S'=Db'JCb so S=Cb'J'Db
%B=B
%R=Db'JDb
%Q=Cb'JCb
%dare(A,B,Q,R,S,E)

B=[B1 B2];
Cb=[C1;zeros(l,n)];
Db=[D11 D12;eye(l) zeros(l,m)];
J=[eye(p) zeros(p,l); zeros(l,p) eye(l)*-gam^2];

[X,L,F]=dare(A,B,Cb'*J*Cb,Db'*J*Db,Cb'*J'*Db);
F=-F;



R=Db'*J*Db+B'*X*B;



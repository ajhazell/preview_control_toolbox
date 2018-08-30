function [X,F,L]=X2d(sys)
% [X,F,L]=X2d(sys)
% Compute solution to the Hinfinity Full Information Riccati equation
% associated with P
%
% X:    DARE solution
% F:    feedback gain 
% L:    closed-loop eigenvalues

[A,B1,B2,C1,C2,D11,D12,D21,D22]=GetSS(sys);
[n,p,q,l,m]=Getsz(sys);

%E'XE=A'XA-(A'XB+S)(B'XB+R)^-1(A'XB+S)'+Q
%

%dare(A,B,Q,R,S,E)

[X,L,F]=dare(A,B2,C1'*C1,D12'*D12,C1'*D12);

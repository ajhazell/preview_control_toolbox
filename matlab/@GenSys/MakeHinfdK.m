function [K,rep]=MakeHinfdK(sys,gam,FI_only,X,Z)
% Generate Hinf-optimal, discrete-time, output feedback controller
% X is the solution to the full information Hinf Riccati Equation
% Z is the solution to the output feedback Hinf Riccati Equation
% Ref: Green and Limebeer, Linear Robust Control
if nargin==3
	X=Xinfd(sys,gam);
	Y=Xinfd(adj(sys),gam);
	Z=Y*inv(eye(size(X))-gam^-2 * X*Y);
end
if ~FI_only
    FI_only=0;
end

rep=0;
[A,B1,B2,C1,C2,D11,D12,D21,D22]=GetSS(sys);
[n,p,q,l,m]=Getsz(sys);


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
L1=L(1:l,:);
R3inv=inv(R3);


Jh=[eye(l) zeros(l,m); zeros(m,l) -gam^2*eye(m)];
[vr3,dr3]=svd(R3);
V12=sqrt(dr3)*vr3;
if norm(V12'*V12-R3)>1e-6
	error('V12 Decomposition has failed')
end
nab=R1-R2'*R3inv*R2;

%[vnab,dnab]=svd(-gam^(-2)*nab);
%V21=sqrt(dnab)*vnab;
%V21=sqrtm(-gam^(-2)*nab); % don't know why this works!!!
V21=chol(-gam^(-2)*nab);
if norm(V21'*V21+gam^(-2)*nab)>1e-8
	V21
	%dnab
	%vnab
    norm(V21'*V21+gam^(-2)*nab)
	error('V21 Decomposition has failed')
end


'Lnab';
Lnab=L1-R2'*R3inv*L2;

'nabinv';
nabinv=inv(nab);
'A';
A;
'B1';
B1;
'At';
At=A-B1*nabinv*Lnab;
Ct1=V12*R3inv*(L2-R2*nabinv*Lnab);
Ct2=C2-D21*nabinv*Lnab;
'Ct';
Ct=[Ct1;Ct2];

Bttil=[B1*inv(V21) zeros(n,m)];
V21inv=inv(V21);
V12inv=inv(V12);
Dttil=[V12*R3inv*R2*V21inv  eye(m);D21*V21inv zeros(q,m)];

St=Dttil*Jh*Dttil'+Ct*Z*Ct';
Mt=Bttil*Jh*Dttil'+At*Z*Ct';

St1=St(1:m,1:m);
St2=St(1:m,1+m:end);
St3=St(1+m:m+q,1+m:m+q);
Mt2=Mt(:,m+1:m+q);
St3inv=inv(St3);



AK=At-Mt2*St3inv*Ct2-B2*V12inv*(Ct1-St2*St3inv*Ct2);
BK=Mt2*St3inv-B2*V12inv*St2*St3inv;
CK=-V12inv*Ct1+V12inv*St2*St3inv*Ct2;
DK=-V12inv*St2*St3inv;

%K=ss(AK,BK,CK,DK,sys.ss.Ts);

if FI_only
    K=-inv(R3)*[L2 R2];
else
    CKt=inv(eye(m)+DK*D22)*CK;
    DKt=inv(eye(m)+DK*D22)*DK;
    AKt=AK-BK*D22*CKt;
    BKt=BK-BK*D22*DKt;
    K=ss(AKt,BKt,CKt,DKt,sys.ss.Ts);
end

% 
% Pm=pck(A,[B1 B2],[C1;C2],[D11 D12; D21 D22]);
% [Km,CLm,gamm,info] = dhfsyn(Pm,q,m,gam,gam,0.001,Ts);
% size(Km)
% size(K)
% 'check K is right ( this should be small)'
% norm(Km-K,inf)

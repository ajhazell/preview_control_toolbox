function K=MakeH2dK(P,X,Y)
% K=MakeH2dK(P,X,Y)
%
% Generate H2-optimal, discrete-time, output feedback controller
% X is the solution to the full information H2 Riccati Equation
% Y is the solution to the output feedbakc H2 Riccati Equation

[A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(P);

Rb=B2'*X*B2+D12'*D12;
F2=-inv(Rb)*(B2'*X*A+D12'*C1);
F0=-inv(Rb)*(B2'*X*B1+D12'*D11);


T=D21*D21'+C2*Y*C2';
L2=-(A*Y*C2'+B1*D21')*inv(T);
L0=(F2*Y*C2'+F0*D21')*inv(T);

%compute controller assuming D22=0
AK=A+B2*F2+L2*C2-B2*L0*C2;
BK=-L2+B2*L0;
CK=F2-L0*C2;
DK=L0;


%Account for non zero D22
CKt=inv(eye(size(B2,2))+DK*D22)*CK;
DKt=inv(eye(size(B2,2))+DK*D22)*DK;
AKt=AK-BK*D22*CKt;
BKt=BK-BK*D22*DKt;

K=ss(AKt,BKt,CKt,DKt,P.ss.Ts);




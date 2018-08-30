function K=MakeH2dK(P)
% Generate H2-optimal, discrete-time, output feedback controller
% X is the solution to the full information H2 Riccati Equation
% Y is the solution to the output feedbakc H2 Riccati Equation
[Ag,B1gw,B2g,C1g,C2g,D11g,D12g,D21g,D22g]=GetSS(P.GW); % ss matrices for G with Wz and Ww absorbed
[B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(P.GW);
[Ar,Br,Cr,Dr]=ssdata(P.Wr);
[n,p,q,l,m]=Getsz(P);
[ng,p,qg,l,m]=Getsz(P.GW);
lr=P.lr;
lw=P.lw;
nwr=P.nwr;


[K,F2,F0]=KFI2d(P)

Sgg=D21gw*D21gw'+C2g*Ygg*C2g';
L2g=-(Ag*Ygg*C2g'+B1gw*D21gw')*inv(Sgg);
L0y=(F2g*Ygg*C2g'+F0w*D21gw')*inv(Sgg);

%compute controller assuming D22=0
AKgg=Ag+B2g*F2g+L2g*C2g-B2g*L0y*C2g;
BKgy=-L2g+B2g*L0y;
CKg=F2g-L0y*C2g;
DKy=L0y;



%Account for non zero D22
CKt=inv(eye(size(B2,2))+DK*D22)*CK;
DKt=inv(eye(size(B2,2))+DK*D22)*DK;
AKt=AK-BK*D22*CKt;
BKt=BK-BK*D22*DKt;

K=ss(AKt,BKt,CKt,DKt,P.ss.Ts);



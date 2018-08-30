function [AKgg,AKgp,AKgr,BKgy,BKgr,CKg,CKp,CKr,L0y,L0r]=GetH2dKParam(P)

[Ag,B1gw,B2g,C1g,C2g,D11g,D12g,D21g,D22g]=GetSS(P.GW); % ss matrices for G with Wz and Ww absorbed
[B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(P.GW);
[Ar,Br,Cr,Dr]=ssdata(P.Wr);
[n,p,q,l,m]=Getsz(P);
[ng,p,qg,l,m]=Getsz(P.GW);
[Ar,Br,Cr,Dr]=ssdata(P.Wr);
lr=P.lr;
lw=P.lw;
nwr=P.nwr;
nps=P.N*P.lr;
Cp=[eye(lr) zeros(lr,nps-lr)];

[K,F2,F0]=KFI2d(P);

F2g=F2(:,1:ng);
F2p=F2(:,ng+1:ng+nps);
F2r=F2(:,ng+nps+1:end);
F0r=F0(:,1:lr);
F0w=F0(:,lr+1:end);

Ygg=X2d(adj(P.GW(:,lr+1:end)));

Sgg=D21gw*D21gw'+C2g*Ygg*C2g';
L2g=-(Ag*Ygg*C2g'+B1gw*D21gw')*inv(Sgg);
L0y=(F2g*Ygg*C2g'+F0w*D21gw')*inv(Sgg);
L0r=F0r*inv(Dr);

%compute controller assuming D22=0
AKgg=Ag+B2g*F2g+L2g*C2g-B2g*L0y*C2g;
AKgp=B1gr*Cp+B2g*F2p+L2g*D21gr*Cp-B2g*L0y*D21gr*Cp;
AKgr=B2g*F2r-B2g*F0r*inv(Dr)*Cr;
BKgy=-L2g+B2g*L0y;
BKgr=B2g*F0r*inv(Dr);
CKg=F2g-L0y*C2g;
CKp=F2p-L0y*D21gr*Cp;
CKr=F2r-F0r*inv(Dr)*Cr;
DKy=L0y;
DKr=L0r;




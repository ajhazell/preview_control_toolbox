function K=MkH2PrevCont(P);
% K=MkH2PrevCont(P);
%
% Construct preview controller
% 

[A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(P);
%[X,F,L]=X2d(P);
%Y=X2d(adj(P));
[AKgg,AKgp,AKgr,BKgy,BKgr,CKg,CKp,CKr,L0y,L0r]=GetH2dKParam(P);
nps=P.N*P.lr;

Ap=[zeros(nps-P.lr,P.lr) eye(nps-P.lr,nps-P.lr);zeros(P.lr,nps)];
Bp=[zeros(nps-P.lr,P.lr);eye(P.lr)];
Cp=[eye(P.lr) zeros(P.lr,nps-P.lr)];

AK=[AKgg AKgp AKgr;zeros(nps,P.GW.n) Ap zeros(nps,P.nwr);zeros(P.nwr,nps+P.GW.n), P.Wr.A-P.Wr.B*inv(P.Wr.D)*P.Wr.C];
BK=[BKgy BKgr;zeros(nps,P.GW.q) Bp; zeros(P.nwr,P.GW.q) P.Wr.B*inv(P.Wr.D)];
CK=[CKg CKp CKr];
DK=[L0y,L0r];

%Account for non zero D22
CKt=inv(eye(size(B2,2))+DK*D22)*CK;
DKt=inv(eye(size(B2,2))+DK*D22)*DK;
AKt=AK-BK*D22*CKt;
BKt=BK-BK*D22*DKt;
Ptmp=P;
K=ss(AKt,BKt,CKt,DKt,P.DistRejGSys.Ts);
%K=PrevController(K,P.lr,P.N);



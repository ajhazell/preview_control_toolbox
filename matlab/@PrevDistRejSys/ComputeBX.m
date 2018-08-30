function [BX]=ComputeBX(P,N,gam)
% BX=ComputeBX(P,N,gam)
%
% Computes B'X. Used as a helper function for MkHinfPrevK

rep=0;

if N<1
	error('PrevTools:ntoosmall','Need N>0')
end
[ng,p,qg,l,m]=Getsz(P.GW);
lw=P.lw;
lr=P.lr;
nwr=P.nwr;
nps=N*lr;
[Ar,Br,Cr,Dr]=ssdata(P.Wr);
[A1prev,B1prev,R1prev,L1prev,Q1prev]=GetSSb(SetN(P,1),gam); 
hs=1:ng+lr;
gs=[1+lr:lr+lw+m];
Ahk=A1prev(hs,hs);
Bhk=[B1prev(1:ng,gs); zeros(lr,m+lw)];
Qhk=Q1prev(hs,hs);
Lhk=L1prev(hs,gs);
Rhk=R1prev(gs,gs);

Ar=A1prev(end-nwr+1:end,end-nwr+1:end);
BpCr=A1prev(1:end-nwr,end-nwr+1:end);
BtilNm1=B1prev(end-nwr-lr+1:end,:);
tic

Lambdak=zeros(lw+m,ng);
betahk=eye(ng);
Psibgk=zeros(ng);%DEBUGGING

for k=1:N-1
	Bhkm1=Bhk;
	[Ahk,Bhk,Qhk,Rhk,Lhk,Hhk,T1hk]=Reduce(P,Ahk,Bhk,Qhk,Rhk,Lhk,ng,lr,nwr);
	Hh11k=Hhk(1:ng,1:ng);
    Hh12k=Hhk(1:ng,1+ng:ng+lr);
    Lambdak=[Lambdak+Bhkm1(1:ng,:)'*Hh11k*betahk', Bhkm1(1:ng,:)'*Hh12k];
	Psibgk=[Psibgk+betahk(1:ng,:)*Hh11k*betahk', betahk(1:ng,:)*Hh12k];%DEBUGGING
	T11hk=T1hk(1:ng,:);
    T12hk=T1hk(ng+1:ng+lr,:);
	betahk=[betahk*T11hk;T12hk];
	
end

ANm1=[Ahk BpCr;zeros(nwr,ng+lr) Ar];
BNm1=[zeros(ng,lr)  Bhk(1:ng,:);BtilNm1];
QNm1=[Qhk zeros(ng+lr,nwr);zeros(nwr,nwr+ng+lr)];
RNm1=[ -eye(lr)*gam^2 zeros(lr,m+lw);zeros(m+lw,lr) Rhk];
LNm1=[zeros(ng+lr,lr) Lhk;zeros(nwr,m+lw+lr)];
%LambdaNm1=[Lambdak zeros(lw+m,nwr+lr)];

%betaNm1=[betahk zeros(ng+nps-lr,nwr+lr);zeros(nwr+lr,ng) eye(nwr+lr)];
[AN,BN,QN,RN,LN,HNm1,T1Nm1]=ReduceFull(P,ANm1,BNm1,QNm1,RNm1,LNm1,ng,lr,nwr);
Hh11Nm1=HNm1(1:ng,1:ng);
Hh12Nm1=HNm1(1:ng,1+ng:ng+lr);
Hh22Nm1=HNm1(1+ng:ng+lr,1+ng:ng+lr);
LambdaN=[Lambdak+Bhk(1:ng,:)'*Hh11Nm1*betahk', Bhk(1:ng,:)'*Hh12Nm1];

T11hNm1=T1Nm1(1:ng,1:ng);
T12hNm1=T1Nm1(ng+1:ng+lr,1:ng);
betahN=[betahk*T11hNm1;T12hNm1];
betaN=[betahN zeros(ng+nps,nwr);zeros(nwr,ng) eye(nwr)];
try
	XN=dare(AN,BN,QN,RN,LN);
catch
	error('PrevTools:XNm1compfailed','Failed to compute XNm1 during reduction algorithm')
end

BX=[Dr'*Hh12Nm1'*betahk' , Dr'*Hh22Nm1, zeros(lr,nwr);LambdaN zeros(lw+m,nwr)]+BN'*XN*betaN';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5




function [Xg,Xr,Rb]=ComputeXgXr(P,N,gam)
% [Xg,Xr,Rb]=ComputeXgXr(P,N,gam)
%
% Computes Xg and Xr (notation from PhD Thesis)
% Acts as helper function for MkHinfPrevK

rep=0;

if N<1
	error('PrevTools:ntoosmall','Need N>0')
end
[ng,p,qg,l,m]=Getsz(P.GW);
lw=P.lw;
lr=P.lr;
nwr=P.nwr;
nps=N*lr;

[A1prev,B1prev,R1prev,L1prev,Q1prev]=GetSSb(SetN(P,1),gam); 
hs=1:ng+lr;
gs=[1+lr:lr+lw+m];
Ahk=A1prev(hs,hs);
Bhgk=[B1prev(1:ng,gs); zeros(lr,m+lw)];
Qhk=Q1prev(hs,hs);
Lhgk=L1prev(hs,gs);
Rgk=R1prev(gs,gs);

Ar=A1prev(end-nwr+1:end,end-nwr+1:end);
BpCr=A1prev(1:end-nwr,end-nwr+1:end);
BtilNm1=B1prev(end-nwr-lr+1:end,:);
tic

Psibgk=zeros(ng);
betabk=eye(ng);
for k=1:N-1
	
	[Ahk,Bhgk,Qhk,Rgk,Lhgk,Hhk,T1bk]=Reduce(Ahk,Bhgk,Qhk,Rgk,Lhgk,ng,lr,nwr);
	Hh11k=Hhk(1:ng,1:ng);
    Hh12k=Hhk(1:ng,1+ng:ng+lr);
    Psibgk=[Psibgk+betabk(1:ng,:)*Hh11k*betabk', betabk(1:ng,:)*Hh12k];
	T11bk=T1bk(1:ng,:);
    T12bk=T1bk(ng+1:ng+lr,:);
	betabk=[betabk*T11bk;T12bk];
	
end

ANm1=[Ahk BpCr;zeros(nwr,ng+lr) Ar];
BNm1=[zeros(ng,lr)  Bhgk(1:ng,:);BtilNm1];
QNm1=[Qhk zeros(ng+lr,nwr);zeros(nwr,nwr+ng+lr)];
RNm1=[ -eye(lr)*gam^2 zeros(lr,m+lw);zeros(m+lw,lr) Rgk];
LNm1=[zeros(ng+lr,lr) Lhgk;zeros(nwr,m+lw+lr)];
PsigNm1=[Psibgk zeros(ng,nwr+lr)];

betaNm1=[betabk zeros(ng+nps-lr,nwr+lr);zeros(nwr+lr,ng) eye(nwr+lr)];

try
	XNm1=dare(ANm1,BNm1,QNm1,RNm1,LNm1);
catch
	error('PrevTools:XNm1compfailed','Failed to compute XNm1 during reduction algorithm')
end

Xg=PsigNm1+betaNm1(1:ng,:)*XNm1*betaNm1';
Xr=betaNm1(end-lr-nwr+1:end,:)*XNm1*betaNm1';
Rb=RNm1+BNm1'*XNm1*BNm1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5


function [Ahkp1,Bhgkp1,Qhkp1,Rgkp1,Lhgkp1,Hhk,T1b,rep]=Reduce(Ahk,Bhgk,Qhk,Rgk,Lhgk,ng,lr,nwr)
rep=0;

U=null([Ahk Bhgk;Lhgk' Rgk]);



try
	U=U(:,1:lr); % just use lr cols of null space
catch
	error('PrevTools:nullspacecompfailed','Failed to compute null space during reduction algorithm')
end

[u,s,v]=svd(U(1:ng+lr,:),0);
P=v*inv(s);
T2b=u;
T1b=null(T2b');
Tb=[T1b T2b];

T1h=[T1b zeros(ng+lr,lr)];
T2h=[T2b zeros(ng+lr,lr)];

Vb=Qhk*T2b+Lhgk*U(1+ng+lr:end,:)*P;
Hhk=T2b*Vb'*T1b*T1b'+Vb*T2b';

Abk=[Ahk*T1b [zeros(ng,lr);eye(lr)]];
Qhkp1=T1h'*Qhk*T1h+Abk'*Hhk*Abk;
Rgkp1=Rgk+Bhgk'*Hhk*Bhgk;
Lhgkp1=T1h'*Lhgk+Abk'*Hhk*Bhgk;
Ahkp1=T1h'*Abk;
Bhgkp1=T1h'*Bhgk;





function [Xgg,Rb]=ComputeXggRb(P,N,gam)
% [Xgg,Rb]=ComputeXggRb(P,N,gam)
%
% Computes Xgg and Rb (notation from PhD Thesis).
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

Psiggk=zeros(ng);
betagk=eye(ng);
for k=1:N-1
	
	[Ahk,Bhgk,Qhk,Rgk,Lhgk,Hhk,T1bk]=Reduce(P,Ahk,Bhgk,Qhk,Rgk,Lhgk,ng,lr,nwr);
	Hh11k=Hhk(1:ng,1:ng);

    Psiggk=Psiggk+betagk*Hh11k*betagk';
	T11bk=T1bk(1:ng,:);

	betagk=betagk*T11bk;
	
end


ANm1=[Ahk BpCr;zeros(nwr,ng+lr) Ar];
BNm1=[zeros(ng,lr)  Bhgk(1:ng,:);BtilNm1];
QNm1=[Qhk zeros(ng+lr,nwr);zeros(nwr,nwr+ng+lr)];
RNm1=[ -eye(lr)*gam^2 zeros(lr,m+lw);zeros(m+lw,lr) Rgk];
LNm1=[zeros(ng+lr,lr) Lhgk;zeros(nwr,m+lw+lr)];

[AN,BN,QN,RN,LN,HN,T1Nm1]=ReduceFull(P,ANm1,BNm1,QNm1,RNm1,LNm1,ng,lr,nwr);
Hh11N=HN(1:ng,1:ng);
T11bk=T1Nm1(1:ng,1:ng);   
Psiggk=Psiggk+betagk*Hh11N*betagk';

betagN=betagk*T11bk;
PsiggN=Psiggk;


try
	XN=dare(AN,BN,QN,RN,LN);
catch
	error('PrevTools:XNm1compfailed','Failed to compute XNm1 during reduction algorithm')
end

Xgg=PsiggN+betagN*XN(1:ng,1:ng)*betagN';
Rb=RN+BN'*XN*BN ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5





%*************************************************************
function [Akp1,Bkp1,Qkp1,Rkp1,Lkp1,Hk,T1]=ReduceFull(P,Ak,Bk,Qk,Rk,Lk,ng,lr,nwr)
nk=size(Ak,1);
hs=1:ng+lr;
%gs=[lr+1:end];
Ahk=Ak(hs,hs);
Bhgk=[Bk(1:ng+lr,lr+1:end)];
Qhk=Qk(hs,hs);
Lhgk=Lk(hs,lr+1:end);
Rgk=Rk(lr+1:end,lr+1:end);
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

T1=[T1b,zeros(ng+lr,nk-lr-ng);zeros(nk-ng-lr,ng) eye(nk-lr-ng)];
T2=[T2b;zeros(nk-ng-lr,lr)];
W=[zeros(lr,lr) ;U(1+ng+lr:end,:)*P];
V=Qk*T2+Lk*W;
Hk=T2*V'*T1*T1'+V*T2';


Qkp1=T1'*Qk*T1+T1'*Ak'*Hk*Ak*T1;
Rkp1=Rk+Bk'*Hk*Bk;
Lkp1=T1'*Lk+T1'*Ak'*Hk*Bk;
Akp1=T1'*Ak*T1;
Bkp1=T1'*Bk;
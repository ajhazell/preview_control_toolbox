function [Ahkp1,Bhgkp1,Qhkp1,Rgkp1,Lhgkp1,Hhk,T1b,rep]=Reduce(P,Ahk,Bhgk,Qhk,Rgk,Lhgk,ng,lr,nwr)
rep=0;

U=null([Ahk Bhgk;Lhgk' Rgk]);



try
	U=U(:,1:lr); % just use lr cols of null space
catch
	error('PrevTools:nullspacecompfailed','Failed to compute null space during reduction algorithm')
end

[u,s,v]=svd(U(1:ng+lr,:),0);
P=v*inv(s);% orthogonalising matrix
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

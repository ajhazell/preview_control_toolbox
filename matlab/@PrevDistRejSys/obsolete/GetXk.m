function [Xk,rep]=GetXk(P,N,gamma)
% [Xk,rep]=GetXk(P,N,gamma)
% Efficient method for computing solution to the Full Information
% Hinfinity Riccati equation 
% **EXPERIMENTAL**
tic;
[n,p,q,l,m]=Getsz(Set_k(P,1));
[lw,lr]=Getlwlr(P);
ng=n-lr;
n1=n;
clear n;

[A_1,B1_1,B2_1,C1_1,C2_1,D11_1,D12_1,D21_1,D22_1,Ts]=GetSS(Set_k(P,1));

k=1;

B_1=[B1_1 B2_1];
Bk=B_1;

E_1=[zeros(ng,lr);eye(lr)];
Ek=E_1;
Eh_1=[eye(ng);zeros(lr,ng)];
Ehk=Eh_1;
Db=[D11_1 D12_1;eye(l) zeros(l,m)];
Cb_1=[C1_1;zeros(l,ng+lr)];

J=[eye(p) zeros(p,l); zeros(l,p) eye(l)*-gamma^2];

[Xk,t1,t2,repx]=dare(A_1,B_1,Cb_1'*J*Cb_1,Db'*J*Db,Cb_1'*J'*Db);
TXk=Xk;
if repx<0
    X=0;
    rep=-1;
    disp(['Failed to compute X1 for gamma=' num2str(gamma)])
    return;
end    
%  F_1=-inv(Db'*J*Db+B_1'*Xk*B_1)*((Cb_1'*J*Db)'+B_1'*Xk*A_1);
%  Ack=A_1+B_1*F_1;

for k=1:N-1

Pk=Ek-Bk*inv(Db'*J*Db+B_1'*TXk(:,[1:ng  end-lr+1:end])*B_1)*B_1'*(TXk(:,[1:ng  end-lr+1:end]))*E_1;    
Wk=gamma^-2*E_1'*TXk*Pk+eye(lr);
if min(abs(eig(Wk)))<1e-6
    X=0;
    rep=-2;
    disp('Wk is poorly conditioned')
    return
end


%Xk=[Xk-Ack'*Xk*Ek*inv(Wk)*Ek'*Ack  Ack'*Xk*Ek*inv(Wk);(Ack'*Xk*Ek*inv(Wk))' Ek'*Xk*Pk*inv(Wk)];
% X1k=[eye(ng+k*lr) zeros(ng+k*lr,lr);Ek'*Ack Wk];
% X1kinv=[eye(ng+k*lr) zeros(ng+k*lr,lr);-inv(Wk)*Ek'*Ack inv(Wk)];

Fk=-inv(Db'*J*Db+B_1'*TXk(:,[1:ng  end-lr+1:end])*B_1)*(([Cb_1 zeros(size(Cb_1,1),k*lr-lr)]'*J*Db)'+...
	B_1'*[TXk(:,1:ng+lr)*A_1(1:ng+lr,1:ng+lr) TXk(:,ng+1:end-lr)]);
EhtAck=[A_1(1:ng,:) zeros(ng,k*lr-lr)]+(Ehk'*Bk)*Fk; % Ehk'*Ack
EtAck=(Ek'*Bk)*Fk;%Ek'*Ack
EtXAck=TXk(end-lr+1:end,1:ng)*EhtAck+TXk(end-lr+1:end,end-lr+1:end)*EtAck+[zeros(lr,ng+lr) TXk(end-lr+1:end,1+ng:end-lr)];
%t1=toc;

TXk11=Eh_1'*TXk-Ehk'*EtXAck'*inv(Wk)*EtAck;
TXk21=inv(Wk)*EtXAck;
TXk12=TXk21(:,1:ng)';
TXk22=E_1'*TXk*Pk*inv(Wk);
TXk=[TXk11 TXk12;TXk21 TXk22];
max(svd(EtXAck'*EtAck)) % This should go to zero as gamma increase
max(svd(EtXAck))
max(svd(EtAck))
%t2=toc;


% Wk
% inv(Wk)*gamma^2
%disp(['min(abs(eig(Wk)))=' num2str(min(abs(eig(Wk))))])
%norm(Testmat,inf)
EtXAcktinvWk=EtXAck'*inv(Wk);
% [Av,H,V]=GetSSv(Set_k(P,k),gamma);
% Xk-EtXAcktinvWk*EtAck
% (-H+Ack'*Xk*Av) 
Xk=[Xk-EtXAcktinvWk*EtAck  EtXAcktinvWk;EtXAcktinvWk' Xk(end-lr+1:end,:)*Pk*inv(Wk)];
%t3=toc;
 %'check Xk ok'
 %norm(Xk-Xinfd(Set_k(P,k+1),gamma),inf)

% Ack=X1k*[Ack Pk;zeros(lr,ng+k*lr+lr)]*X1kinv;

Ek=[zeros(lr,lr);Ek];
Ehk=[Ehk;zeros(lr,ng)];

Bk=[B_1(1:end-lr,:);zeros(k*lr,l+m);B_1(end-lr+1:end,:)];

end
%disp(['Low storage time = ' num2str(t2-t1)])
%disp(['High storage time = ' num2str(t3-t2)])
rep=0;
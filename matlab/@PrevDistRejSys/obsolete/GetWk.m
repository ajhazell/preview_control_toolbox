function Wk=GetWk(P,gamma)
% Get kth W - a parameter used in recursive construction of the 
% deflating subspace associated with Extended Symplectic Pencil
% corresponding to the Hinfinty Riccati equation
[n,p,q,l,m]=Getsz(P);
[ng,pg,qg,lg,mg]=Getsz(GetWG(P));
[lw,lr]=Getlwlr(P);

%[A_1,B1_1,B2_1,C1_1,C2_1,D11_1,D12_1,D21_1,D22_1,Ts]=GetSS(Set_k(P,1));
[Ak,B1k,B2k,C1k,C2k,D11k,D12k,D21k,D22k,Ts]=GetSS(P);
 
%B_1=[B1_1 B2_1];
Bk=[B1k B2k];

%E_1=B1_1(:,lw+1:lr+lw);
Ek=B1k(:,lw+1:lr+lw);

Db=[D11k D12k;eye(l) zeros(l,m)];
J=[eye(p) zeros(p,l); zeros(l,p) eye(l)*-gamma^2];

%TXk=GetXk(P,P.k,gamma);
Xk=Xinfd(P,gamma);
%Pk=Ek-Bk*inv(Db'*J*Db+B_1'*TXk(:,[1:ng  end-lr+1:end])*B_1)*B_1'*(TXk(:,[1:ng  end-lr+1:end]))*E_1;    
Pk=Ek-Bk*inv(Db'*J*Db+Bk'*Xk*Bk)*Bk'*Xk*Ek;    
Wk=gamma^-2*Ek'*Xk*Pk+eye(lr);

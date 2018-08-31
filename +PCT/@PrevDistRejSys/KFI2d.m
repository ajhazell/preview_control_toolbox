function [K,F2,F0]=KFI2d(P,N)
% [K,F2,F0]=KFI2d(P)
%
% Efficiently generate full information, discrete-time, H2-optimal  controller gain
% Value of N overrides value stored in P.N
if nargin==1
    N=P.N;
end

[Ag,B1gw,B2g,C1g,C2g,D11g,D12g,D21g,D22g]=GetSS(P.GW); % ss matrices for G with Wz and Ww absorbed
[B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(P.GW);
[Ar,Br,Cr,Dr]=ssdata(P.Wr);
[n,p,q,l,m]=Getsz(P);
[ng,p,qg,l,m]=Getsz(P.GW);
lr=P.lr;
lw=P.lw;
nwr=P.nwr;
D12=D12g;

Xgg=X2d(P.GW);


Rb=B2g'*Xgg*B2g+D12'*D12;
F2g=-inv(Rb)*(B2g'*Xgg*Ag+D12'*C1g);
F0w=-inv(Rb)*(B2g'*Xgg*B1gw+D12'*D11gw);

S=Ag'*Xgg*B1gr+F2g'*B2g'*Xgg*B1gr+F2g'*D12'*D11gr+C1g'*D11gr;
Acg=Ag+B2g*F2g;

Xgp=zeros(ng,(P.N)*P.lr);
F2p0=-inv(Rb)*(B2g'*Xgg*B1gr+D12g'*D11gr);


	
Xgp(:,1:lr)=S;
for i=1:N-1
    Xgp(:,i*lr+[1:lr])=Acg'*Xgp(:,(i-1)*lr+[1:lr]);
end
F2p=[F2p0,-inv(Rb)*B2g'*Xgp(:,1:(N-1)*lr)];

Xgr=dlyap(Acg',Ar,(Acg')*Xgp(:,(end-lr+1):end)*Cr);


F2r=-inv(Rb)*(B2g'*Xgp(:,(end-lr+1):end)*Cr+B2g'*Xgr*Ar);

F0r=-inv(Rb)*(B2g'*Xgr*Br+B2g'*Xgp(:,(end-lr+1):end)*Dr);
F2=[F2g F2p F2r];
F0=[F0r F0w];
K=[F2 F0];
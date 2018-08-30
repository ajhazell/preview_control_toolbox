function [gam,gamprev]=ComputeMinFICL2norm(P)

%[gam,gamprev]=ComputeMinFICL2norm(P)
%
% Computes the minimum ahchievable H2 norm for the plant P
% with preview length P.N.
%
% gam:      min FI 2-norm
% gamprev:  reduction in 2-norm which results from using preview

[Ag,B1gw,B2g,C1g,C2g,D11g,D12g,D21g,D22g]=GetSS(P.GW); % ss matrices for G with Wz and Ww absorbed
[B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(P.GW);
[Ar,Br,Cr,Dr]=ssdata(P.Wr);
[n,p,q,l,m]=Getsz(P);
[ng,p,qg,l,m]=Getsz(P.GW);
lr=P.lr;
lw=P.lw;
nwr=P.nwr;
D12=D12g;
N=P.N;
Xgg=X2d(P.GW);

if P.nwr>0
    error('PrevTools:BadP','Wr is not yet supported for efficient norm computation')
end

Rb=B2g'*Xgg*B2g+D12'*D12;
F2g=-inv(Rb)*(B2g'*Xgg*Ag+D12'*C1g);
F0w=-inv(Rb)*(B2g'*Xgg*B1gw+D12'*D11gw);
F2p0=-inv(Rb)*(B2g'*Xgg*B1gr+D12g'*D11gr);

S=Ag'*Xgg*B1gr+F2g'*B2g'*Xgg*B1gr+F2g'*D12'*D11gr+C1g'*D11gr;
Acg=Ag+B2g*F2g;

[K,F2,F0]=KFI2d(P);

cln=zeros(lr);
for i=0:N-1
    cln=cln-S'*Acg^i*B2g*inv(Rb)*B2g'*(Acg^i)'*S;
end
gamprev=sqrt(trace(-cln));
cln=cln+(B1gr'*Xgg*B1gr+D11gr'*D11gr-F2p0'*Rb*F2p0);
gam=sqrt(trace(cln));


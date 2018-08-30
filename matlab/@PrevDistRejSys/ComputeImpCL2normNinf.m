function [gamprevNinfsq]=ComputeImpCL2normNinf(P)
%[gamprevNinf]=ComputeImpCL2normNinf(P)
%
% Computes the maximum improvement in the squared H2 norm due to preview for the plant P
% with infinite preview length.
%
% gamprevNinfsq:  reduction in squared 2-norm which results from using infinite preview


[Ag,B1gw,B2g,C1g,C2g,D11g,D12g,D21g,D22g]=GetSS(P.GW); % ss matrices for G with Wz and Ww absorbed
[B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(P.GW);
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

S=Ag'*Xgg*B1gr+F2g'*B2g'*Xgg*B1gr+F2g'*D12'*D11gr+C1g'*D11gr;
Acg=Ag+B2g*F2g;

Gamma=dlyap(Acg,B2g*inv(Rb)*B2g');
gamprevNinfsq=(trace(S'*Gamma*S));

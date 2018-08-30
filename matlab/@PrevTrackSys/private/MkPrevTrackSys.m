function P=MkPrevTrackSys(psys)
% Private method used to generate the augmented state space representation 
% of a preview tracking system
[n,p,q,lw,m]=Getsz(psys.G);
lr=psys.lr;
N=psys.N;
l=lw+lr;

[As,B1s,B2s,C1s,C2s,D11s,D12,D21s,D22s,Ts]=GetSS(psys.G);

nps=N*lr;

Ap=[zeros(nps-lr,lr) eye(nps-lr,nps-lr);zeros(lr,nps)];
Bp=[zeros(nps-lr,lr);eye(lr)];

C1=[C1s [ -eye(lr);zeros(p-lr,lr)] zeros(p,nps-lr)];

C2=[C2s zeros(q,nps);zeros(lr,n+nps)];

D21=[D21s zeros(q,lr);zeros(lr,lw) eye(lr)];

D11=[D11s zeros(p,lr)];

D22=[D22s; zeros(lr,m)];

B=[B1s zeros(n,lr) B2s;zeros(nps,lw) Bp zeros(nps,m)];

A=[As zeros(n,nps);zeros(nps,n) Ap];
C=[C1;C2];
D=[D11 D12; D21 D22];

Wwru=ss(eye(l+m));
Wwru(1:lw,1:lw)=psys.Ww;
Wwru(lw+1:l,lw+1:l)=psys.Wr;
Wyz=ss(eye(p+q+lr));
Wyz(1:p,1:p)=psys.Wz;

if strcmp(computer, 'PCWIN')
    P=GenSys(Wyz*ss(A,B,C,D,psys.G.Ts)*Wwru,q+lr,m);
else
    P=GenSys(Wyz*ss(A,B,C,D,psys.G.Ts)*Wwru,q+lr,m);
end


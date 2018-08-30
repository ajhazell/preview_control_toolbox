function [P,psys]=MkPrevDistRejSys(psys)
% Private method used to of generate the augmented state space representation 
% of a preview system
[ng,p,q,l,m]=Getsz(psys.G);
lr=psys.G.lr;
lw=psys.G.lw;
N=psys.N;

% absorb Wz and Ww into G to give GW
Wwru=ss(eye(l+m));
Wwru(lr+1:l,1+lr:l)=psys.Ww;
Wyz=ss(eye(p+q));
Wyz(1:p,1:p)=psys.Wz;
GW=DistRejGSys(Wyz*psys.G*Wwru,q,m,lr);
psys.GW=GW;
[Ag,B1g,B2g,C1g,C2g,D11g,D12,D21g,D22g,Ts]=GetSS(GW);
[B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(GW);

ng=size(Ag,1); % redefine ng

nps=N*lr;

Ap=[zeros(nps-lr,lr) eye(nps-lr,nps-lr);zeros(lr,nps)];
Bp=[zeros(nps-lr,lr);eye(lr)];
Cp=[eye(lr) zeros(lr,nps-lr)];

C1=[C1g D11gr*Cp ];

C2=[C2g D21gr*Cp;zeros(lr,ng+nps)];

D21=[ zeros(q,lr) D21gw;  eye(lr) zeros(lr,lw) ];

D11=[ zeros(p,lr) D11gw];

D22=[D22g; zeros(lr,m)];

B=[zeros(ng,lr) B1gw B2g;Bp zeros(nps,m+lw)];

A=[Ag B1gr*Cp;zeros(nps,ng) Ap];
C=[C1;C2];
D=[D11 D12; D21 D22];

Wwru=ss(eye(l+m));
Wwru(1:lr,1:lr)=psys.Wr;
P=DistRejGSys(ss(A,B,C,D,psys.G.Ts)*Wwru,q+lr,m,lr);


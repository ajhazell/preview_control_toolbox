function P=MkDistRejGSys(gsys)
% [n,p,q,l,m]=Getsz(gsys);
% lr=gsys.lr;
% 
% l=lw+lr;
% 
% [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(gsys);
% 
% C=[C1;C2];
% D=[D11 D12; D21 D22];
% 
% Wwru=ss(eye(l+m));
% Wwru(1:lw,1:lw)=gsys.Ww;
% Wwru(lw+1:l,lw+1:l)=gsys.Wr;
% Wyz=ss(eye(p+q+lr));
% Wyz(1:p,1:p)=gsys.Wz;
% 
% P=GenSys(Wyz*gsys.G*Wwru,q,m);
% 


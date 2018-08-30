function psys=PrevDistRejSys(G,N,Wz,Wr,Ww)
% P=PrevDistRejSys(G,N,Wz,Wr,Ww)
%
% Constructor for PrevSys
% N : Preview length
% Wz : Weighting function for z (can be an lti sys or simple matrix gain)
% Wr : Weighting function for r (can be an lti sys or simple matrix gain)
% Ww : Weighting function for w (can be an lti sys or simple matrix gain)
% G : Object of class DistRejGSys whose output is to be regulated

lw=G.lw;
lr=G.lr;
p=G.p;
qg=G.q;
l=G.l;
m=G.m;
ng=G.n;

psys=struct('lr',lr,'lw',lw,'N',N,'G',G,'ng',ng,'GW',[]);


psys.Wr=ss(eye(lr));
psys.Wz=ss(eye(p));
psys.Ww=ss(eye(lw));
psys.nwr=0;

if nargin>=3;
	if size(ss(Wz))==p
		psys.Wz=ss(Wz);
	else
		error('Incorrect dimensions for Wz')
	end
	
end
	
if nargin>=4;
	if size(ss(Wr))==lr
		psys.Wr=ss(Wr);
		psys.nwr=size(psys.Wr.A,1);
	else
		error('Incorrect dimensions for Wr')
	end
	
end

if nargin==5;
	if size(ss(Ww))==lw
		psys.Ww=ss(Ww);
	else
		error('Incorrect dimensions for Ww')
	end
	
end

if nargin>5 | nargin<2;
    error('Incorrect number of input arguments')
end


[P,psys]=MkPrevDistRejSys(psys);
psys=class(psys,'PrevDistRejSys',P);
superiorto('DistRejGSys')
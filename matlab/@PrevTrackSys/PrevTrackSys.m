function psys=PrevTrackSys(G,N,lr,Wz,Wr,Ww)
% PrevTrackSys(G,N,lr,Wz,Wr,Ww)
%
% Constructor for the class PrevTrackSys, which represents the generalised plant
% associated with the preview tracking problem. PrevTrackSys inherits from
% DistRejPrevSys.
% N : Preview length
% lr : Dimension of previewable reference to be tracked
% Wz : Weighting function for z (can be an lti sys or simple matrix gain)
% Wr : Weighting function for r (can be an lti sys or simple matrix gain)
% Ww : Weighting function for w (can be an lti sys or simple matrix gain)
% G : Object of class GenSys whose is output is to be tracked. G performs
% the mapping [zg;y]=G[w;u]

%note Z=[z1;z2] and z1 is the tracking error term, computed from
%z1=zg-Phi*r

[ng,p,qg,lw,m]=Getsz(G); 

psys.Wr=ss(eye(lr));
psys.Wz=ss(eye(p));
psys.Ww=ss(eye(lw));
psys.nwr=0;

if nargin>=4;
	if size(ss(Wz))==p
		psys.Wz=ss(Wz);
	else
		error('Incorrect dimensions for Wz')
	end
	
end
	
if nargin>=5;
	if size(Wr)==lr
		psys.Wr=ss(Wr);
		psys.nwr=size(psys.Wr.A,1);
	else
		error('Incorrect dimensions for Wr')
	end
	
end

if nargin==6;
	if size(Ww)==lw
		psys.Ww=ss(Ww);
	else
		error('Incorrect dimensions for Ww')
	end
	
end

if nargin>6 | nargin<3;
    error('Incorrect number of input arguments')
end


[ng,p,qg,lw,m]=Getsz(G); % lw is the size of the disturbance w
l=lw+lr;

[Ag,B1gw,B2g,C1g,C2g,D11gw,D12,D21gw,D22g,Ts]=GetSS(G);

C1=[C1g ];

C2=[C2g];

D21=[ zeros(qg,lr)  D21gw  ];
 
D11gr=[ -eye(lr);zeros(p-lr,lr)];
D11=[D11gr D11gw];

D22=[D22g];

B1gr=zeros(ng,lr);
B=[ B1gr B1gw B2g];

A=[Ag ];
C=[C1;C2];
D=[D11 D12; D21 D22];

Gtrack=DistRejGSys(ss(A,B,C,D,G.Ts),qg,m,lr);

P=PrevDistRejSys(Gtrack,N,psys.Wz,psys.Wr,psys.Ww);
psys=class(psys,'PrevTrackSys',P);
superiorto('PrevDistRejSys')





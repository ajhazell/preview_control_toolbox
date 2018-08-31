function PFI=GetFIPlant(P)
% PFI=GetFIPlant(P)
% Generates "Full Information" plant, so y comprises the full state together with 
% the disturbances i.e. y=[x;w]
[A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(P);
PFI=PCT.GenSys(ss(A,[B1 B2],[C1;eye(P.n);zeros(P.l,P.n)],...
	[D11 D12;[zeros(P.n,P.l);eye(P.l)] zeros(P.n+P.l,P.m) ],Ts),P.n+P.l,P.m);

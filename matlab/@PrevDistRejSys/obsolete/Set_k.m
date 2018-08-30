function psys=SetN(psys,N)
% Set the preview length for the Preview System
psys.N=N;
P=MkPrevSys(psys);
psys.GenSys=P;
%psys=class(psys,'PrevSys',P);
%superiorto('GenSys')
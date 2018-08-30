function psys=SetN(psys,N)
% P=SetN(P,N)
% Set the preview length for the Preview System


%psys.N=N;
%P=MkPrevDistRejSys(psys);
%psys.GenSys=P;
%psys=class(struct(psys),'PrevSys',P);
%superiorto('GenSys')
psys=PrevDistRejSys(psys.G,N,psys.Wz,psys.Wr,psys.Ww);
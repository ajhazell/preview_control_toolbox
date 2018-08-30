function psys=Inc_k(psys)
% psys=Inc_k(psys)
% Increase the preview length of psys
psys.k=psys.k+1;
psys=MkPrevSys(psys)
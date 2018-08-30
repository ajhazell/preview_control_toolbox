function WG=GetWG(P)
% WG=GetWG(P)
% Gets weighted transfer func from [w;u] to [z;y]   (i.e. ignores r)

disp('WARNING: Assuming that Ww=I and Wz=I')
%error('Function needs modification to cope with Ww')

[n,p,q,lw,m]=Getsz(P.G);
[lw, lr]=Getlwlr(P);
WG=P.G;

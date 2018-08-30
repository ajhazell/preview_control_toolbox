function gsys=DistRejGSys(G,q,m,lr)
% DistRejGSys(G,q,m,lr)
%
% Constructor for DistRejGSys - Helper class for PrevDistRejSys
%
% q: : Dimension of the measurement (i.e. length(y))
% m  : Dimension of the control (i.e. length(u))
% lr : Dimension of reference to be tracked
% G  : Object of class ss 



l=size(G.B,2)-m;
lw=l-lr;

if lw<0
    error('PrevTools:lwnegative','Inferred value of lw is negative')
end

gsys=struct('lr',lr,'lw',lw);

P=GenSys(G,q,m);
gsys=class(gsys,'DistRejGSys',P);
superiorto('GenSys')


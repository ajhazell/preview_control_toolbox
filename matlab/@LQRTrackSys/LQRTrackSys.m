function psys=LQRTrackSys(G,R,Q,N)
% P=LQRTrackSys(G,R,Q,N)
%   
% Computes generalised plant for LQR preview tracking formulation
%
% Note that the cost function has the form J=sum [e'Qe+u'Ru]
% i.e. Q is an error weighting, not a state weighting
%
% N is the preview length
%
% WARNING: LQRTrackSys is early version code, so its probably quite easy to
% break

m=size(G,2);
p1=size(G,1);
lr=p1;
Gtrack=[G;eye(m)];
We=chol(Q);
Wu=chol(R);
psys=struct('Q',Q,'R',R);

Wz=[We zeros(p1,m);zeros(m,p1) Wu];
P=PrevTrackSys(GenSys(Gtrack,0,m),N,lr,Wz);

psys=class(psys,'LQRTrackSys',P);
superiorto('PrevTrackSys')





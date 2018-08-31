function [n,p,q,l,m]=Getsz(sys);
% [n,p,q,l,m]=Getsz(sys);
% Get sizes of the generalised plant
% 
%  n number of states
%  p number of controlled variables
%  q number of measurements
%  l number of disturbances
%  m number of controls
n=sys.n;
p=sys.p;
q=sys.q;
l=sys.l;
m=sys.m;
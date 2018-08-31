function [B1r,B1w,D11r,D11w,D21r,D21w]=GetSSrw(sys)
% [B1r,B1w,D11r,D11w,D21r,D21w]=GetSSrw(sys)
% 
% Get state space matrix partitions relating to partitioning of the
% exogenous input

q=sys.q;
m=sys.m;

[A,B,C,D,Ts]=ssdata(sys); % discrete time
    
B1=B(:,1:end-m);
D11=D(1:end-q,1:end-m);
D21=D(end-q+1:end,1:end-m);

rsel=1:sys.lr;
wsel=sys.lr+1:sys.l;

B1r=B1(:,rsel);
B1w=B1(:,wsel);

D11r=D11(:,rsel);
D11w=D11(:,wsel);

D21r=D21(:,rsel);
D21w=D21(:,wsel);


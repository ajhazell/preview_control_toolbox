function [A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(sys)
% Get state space matrices for generalised plant
q=sys.q;
m=sys.m;
if nargout==9
    [A,B,C,D]=ssdata(sys); % continuous time
elseif nargout==10
    [A,B,C,D,Ts]=ssdata(sys); % discrete time
else
    error('Incorrect number of output arguments')
end
    
B1=B(:,1:end-m);
B2=B(:,end-m+1:end);
C1=C(1:end-q,:);
C2=C(end-q+1:end,:);
D11=D(1:end-q,1:end-m);
D12=D(1:end-q,end-m+1:end);
D21=D(end-q+1:end,1:end-m);
D22=D(end-q+1:end,end-m+1:end);



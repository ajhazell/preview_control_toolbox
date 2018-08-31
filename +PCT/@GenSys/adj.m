function adjsys=adj(sys)
% adjsys=adj(P)
% OBSELETE FUNCTION DUE TO POOR NAMING - ONLY INCLUDED FOR BACKWARD
% COMPATIBILITY
% adjsys=ss(A',C',B',D',Ts);
% Note: this is NOT the 'adjoint' of  P

[A,B,C,D]=ssdata(sys);
[n,p,q,l,m]=Getsz(sys);

adjsys=ss(A',C',B',D',sys.ss.Ts);
adjsys=PCT.GenSys(adjsys,m,q);
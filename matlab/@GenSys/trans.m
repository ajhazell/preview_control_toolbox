function transys=trans(sys)
% transys=adj(P)
%
% Computes transpose of system P:
% transys=ss(P.A',P.C',P.B',P.D',Ts);
% transys=GenSys(transys,P.m,P.q);


[A,B,C,D]=ssdata(sys);
[n,p,q,l,m]=Getsz(sys);

transys=ss(A',C',B',D',sys.ss.Ts);
transys=GenSys(transys,m,q);
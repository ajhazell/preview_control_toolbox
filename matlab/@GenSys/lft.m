function CL=lft(P,K)
% CL=lft(P,K)
%
% Computes linear fractional transformation of P with K
%
% P: Object of class GenSys
% K: Object of class ss

CL=lft(P.ss,K,P.m,P.q);

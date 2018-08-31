function K=KFI2d(P,X)
% K=KFI2d(P,X)
% Generate full information, discrete-time, H2-optimal  controller gain
[A,B1,B2,C1,C2,D11,D12,D21,D22,Ts]=GetSS(P);
Rb=B2'*X*B2+D12'*D12;
F2=-inv(Rb)*(B2'*X*A+D12'*C1);
F0=-inv(Rb)*(B2'*X*B1+D12'*D11);

K=[F2 F0];
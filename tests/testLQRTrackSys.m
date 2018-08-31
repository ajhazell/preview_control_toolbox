function tests = testLQRTrackSys
    tests = functiontests(localfunctions);
end

% TODO: check closed loop H2 norm is what it should be
% TODO: break up this test case so it tests one part at a time

function tests = testBasicCheckNoErrors(testCase)
% feeble test case that just checks no errors are thrown
% better than nothing though....

% Define plant to be controlled (G)
% Sample Time
Ts=0.002;
cp=-1+10*j; % an stable pole
G=zpk([],[cp,cp'],[100]); % a zero-pole-gain model
G=c2d(G,Ts,'tustin'); % convert to discrete-time
R=1;
Q=10;
N=1000;
P=PCT.LQRTrackSys(G,R,Q,N);
K2=MkH2dK(P); % Compute OF controller, calls efficient routines in DistRejPrevSys
CL2=lft(P,K2); %  form closed loop using a linear fractional transformation

% Compute improvement in 2-norm due to preview
[gamFI2,gam2prev]=ComputeMinFICL2norm(P);

%%
% Compute max possible improvement
[gam2prevmax]=ComputeImpCL2normNinf(P);

% check that our 1000 steps has given us all the benefit we're going to get
verifyEqual(testCase, gam2prev , gam2prevmax, 'RelTol',1e-2)

end
clear classes

%%
% Define signal dimensions
p=2; % dim(z) 
q=1; % dim(y)
m=1; % dim(u)
lr=1;% dim(r)
l=1; % dim(w)=0 in this example

%%
% Preview length
N=200;

%%
% Sample Time
Ts=0.002;

%%
% Define laplace and z-transform variables
z=tf([1 0],1,Ts);
s=tf([1 0],1);

%%
% Define delay line
Phi=z^-N;

%%
% Weightings (for tuning performance)
trackweight=10;
measnoise=1;

%%
% Define plant to be controlled (G)
fprintf('Plant to be controlled:\n')
cp=1+10*j; % an unstable pole
G=zpk([],[cp,cp'],[100]); % a zero-pole-gain model
G=c2d(G,Ts,'tustin'); % convert to discrete-time

%%
% Add additional output so that the control is included in z
Ggen=[ G;  eye(m) ];

%%
% Add an output disturbance
Ggen=[[20;0] Ggen];

%%
% Add a measurement signal (y) - the measured signal is identical to the signal we are trying to control.
Ggen=[ Ggen;Ggen(1,:)];

%%
%Add some measurement noise
Ggen=[[zeros(q+p-1,q); measnoise] Ggen];

%%
% Convert to class GenSys
Ggen=GenSys(minreal(ss(Ggen)),q,m);

%% 
% In this example the outputs of Ggen are [zg+wout;u;y], where zg is the output that should track a delayed version of r, u is the control and y is the measurement, wout is the output disturbance. Also note that in this case, y=zg+noise and also z=Wz*[zg;u]

%% 
% Define low pass output weighting function
Wz=[trackweight*zpk(0,0.99,1,Ts),0;0, 1];

%%
% Define input weighting function
Wr=eye(lr); % could be any 'lti' object

%% 
% Create P, which is the generalised plant for the preview tracking system
P=PrevTrackSys(Ggen,N,lr,Wz,Wr);

%% 
% Notice that PrevTrackSys inherits from DistRejGSys, DistRejPrevSys, and GenSys. It also inherits form the class PublicProperties, which means that one can read its attributes using 'dot' notation e.g. P.Wr, P.A, P.lr, P.Ts etc.

%%
% Compute Output Feedback controller and closed loop
fprintf('Begin computation of OF controllers....')
fprintf('H2...')
K2=MkH2dK(P); % Compute OF controller, calls efficient routines in DistRejPrevSys
CL2=lft(P,K2); %  form closed loop using a linear fractional transformation

fprintf('Hinf...')
gam0=1;h0=1000;hmin=1e-2;bFI=false;
[Kinf,gam]=MkOptHinfdK(P,P.N,gam0,h0,hmin,bFI); % Compute OF controller, calls efficient routines in DistRejPrevSys
fprintf('achieved gamma=%0.3g...',gam)
CLinf=lft(P,Kinf); % linear fractional transformation
fprintf('done.\n')

%Uncomment these lines to do brute force checking of achieved norm
%fprintf('Checking Hinf OF norm bound satisfied...')
%fprintf('Bound=%0.5g, Achieved=%0.5g\n',gam,norm(CLinf,inf))

%%
% Obtain full information version of generalised plant
PFI=GetFIPlant(P);

%%
% Compute Full information controller and closed loop
fprintf('Begin computation of FI controllers....')
fprintf('H2...')
KFI2=KFI2d(P); % Compute FI H2 gains
CLFI2=lft(PFI,KFI2); % linear fractional transformation

fprintf('Hinf...')
gam0=1;h0=1000;hmin=1e-2;bFI=true;
[KFIinf,gamFI]=MkOptHinfdK(P,P.N,gam0,h0,hmin,bFI); % Compute FI controller gain, calls efficient routines in DistRejPrevSys
fprintf('achieved gamma=%0.3g...',gam)
CLFIinf=lft(PFI,KFIinf); % linear fractional transformation, calls efficient routines in DistRejPrevSys
fprintf('done.\n')

%%
% If the preview length is very large, then storing the full augmented
% system, P, may represent an excessive memory requirement. In this
% instance, it is possible to take a P defined for an arbitrary preview length, and
% override this during the computation of the controller. i.e. ignore P.N.
% For example:
Nlong=1000;
fprintf('Computing %5.0f preview point H2 FI controller....',Nlong)
KFI2_long=KFI2d(P,Nlong);
fprintf('done.\n')
%%
% similar overriding is also possible with MkOptHinfdK(...) and MkHinfdK(...)


%%
% Compute improvement in 2-norm due to preview
[gamFI2,gam2prev]=ComputeMinFICL2norm(P);

%%
% Compute max possible improvement
[gam2prevmax]=ComputeImpCL2normNinf(P);

fprintf('Improvement in 2-norm due to preview=%0.3g, which is %5.2f percent of the maximum possible improvement.\n',gam2prev,100*gam2prev/gam2prevmax)


%Uncomment these lines to do brute-force checking of achieved inf norm
%fprintf('Checking Hinf FI norm bound satisfied...')
%fprintf('Bound=%0.5g, Achieved=%0.5g\n',gam,norm(CLinf,inf))

%%
fprintf('\nComputed all controllers, now beginning simulations \n\t[may be slow ... efficient simulation code not yet developed]....')

%%
% Plot FI step responses
Tf=1;
figure(1);clf
subplot(2,1,1)
step(CLFI2(1,1)/(Wz(1,1)*Wr(1,1))+Phi,Tf)
hold all
step(CLFIinf(1,1)/(Wz(1,1)*Wr(1,1))+Phi,Tf)
step(Phi,Tf)
grid on
legend('Ouput of G with H2 Control', 'Output of G with Hinf Control', 'Desired output') 
title('Full Information Preview Tracking Controllers')

subplot(2,1,2)
step(CLFI2(2,1)/(Wz(2,2)*Wr(1,1)),Tf)
hold all
step(CLFIinf(2,1)/(Wz(2,2)*Wr(1,1)),Tf)
grid on
legend('H2 Control','Hinf Control')

%%
% Plot OF step responses
figure(2);clf
subplot(2,1,1)
step(CL2(1,1)/(Wz(1,1)*Wr(1,1))+Phi,Tf)
hold all
step(CLinf(1,1)/(Wz(1,1)*Wr(1,1))+Phi,Tf)
step(Phi,Tf)
grid on
legend('Ouput of G with H2 Control', 'Output of G with Hinf Control', 'Desired output')
title('Output Feedback Preview Tracking Controllers')

subplot(2,1,2)
step(CL2(2,1)/(Wz(2,2)*Wr(1,1)),Tf)
hold all
step(CLinf(2,1)/(Wz(2,2)*Wr(1,1)),Tf)
grid on

fprintf('... done.\n')
function [K,E]=MkHinfPrevK(P,N,gam,bFI)
% [K,E]=MkHinfPrevK(P,N,gam,bFI)
%
% Generate Hinf-suboptimal, discrete-time, output feedback controller
% Ref: Green and Limebeer, Linear Robust Control and A Hazell, Discrete-Time Optimal
% Preview control, PhD Thesis, 2008 
%
% N:    Preview length, overrides the preview length associated with P
% gam:  Maximum allowable hinf norm
% bFI:  boolean, true if only the Full Information (FI) controller is required
% E:    Entropy associated with FI controller
% K:    Optimal FI or OF controller

%reps
%====
%0   :  Success!
%-1  :  Could not compute stabilising X
%-2  :  X was not feasible
%-3  :  X was no positive

rep=0;
tol=1e-9;

if N<1
	error('PrevTools:ntoosmall','Need N>0')
end

if isempty(bFI)
	bFI=false;
end

%%%%%%%%%% Extract partitions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,B1,B2,C1,C2,D11,D12,D21,D22]=GetSS(P);
[Ag,B1gw,B2g,C1g,C2g,D11g,D12g,D21g,D22g]=GetSS(P.GW); % ss matrices for G with Wz and Ww absorbed
[B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(P.GW);
[Ar,Br,Cr,Dr]=ssdata(P.Wr);
[n,p,q,l,m]=Getsz(P);
[ng,p,qg,l,m]=Getsz(P.GW);
lr=P.lr;
lw=P.lw;
nwr=P.nwr;
Lh=[D11';D12']*[C1g D11gr];
Lg=Lh(:,1:ng);


%%%%%%%%%%%%%%%%%%%%Compute Y%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~bFI
	try
		GWg=P.GW(:,lr+1:end);
		GWg=GenSys(ss(GWg.A,GWg.B,GWg.C,GWg.D,GWg.Ts),qg,m);
		Yg=Xinfd(adj(GWg),gam);
		%Yg=0; disp('Warnining: Forcing Yg=0')
		
	catch
		error('PrevTools:Ygcompfailed','Failed to compute Y')
	end

	if max(eig(GetNabla(adj(P.GW(:,lr+1:end)),gam,Yg)))>tol
		error('PrevTools:Ynotfeasible','Y is not feasible')
	end

	if min(eig(Yg))<-tol
		error('PrevTools:Ynotnonnegative','Y is not nonnegative')
	end

end

%%%%%%%%%%%%%%%%%Compute Xgg Rb, and check KFI exists%%%%%%%%%%%%%%%
try
	%[Xg,Xr,Rb]=ComputeXgXr(P,N,gam);
    [Xgg,Rb]=ComputeXggRb(P,N,gam);
catch
	rethrow(lasterror)
end

%Xgg=Xg(:,1:ng);

Rb1=Rb(1:l,1:l);
Rb2=Rb(1:l,l+1:l+m)';
Rb3=Rb(l+1:l+m,l+1:l+m);
Rb3inv=inv(Rb3);

nab=Rb1-Rb2'*Rb3inv*Rb2;
V21=chol(-gam^(-2)*nab);
if norm(V21'*V21+gam^(-2)*nab)>tol
	error('PrevTools:nablanotnegative','Nabla was not negative')
end

Ac_cl_gg=Ag-B2g*Rb3inv*(B2g'*Xgg*Ag+Lg(end-m+1:end,:));
if max(abs(eig(Ac_cl_gg)))>=1
    error('PrevTools:Xnotpositive','X was not non-negative')
end

%%%%%%%%%%%%%%%%%%%%%%%Compute Z %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~bFI
    
    if max(eig(Yg*Xgg))>tol+gam^2
		error('PrevTools:XYconditionnotmet','Spectral radius of XY is greater than gam^2')
	end
    
    Zg=Yg*inv(eye(ng)-gam^-2 * Xgg*Yg);
end

%%%%%%%%%%%%%%%%%%%% Compute B'X %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%If we've got this far then all conditions for existence of K have been
%passed, therefore it is worth computing B'X.

try
	BX=ComputeBX(P,N,gam);
catch
	error('PrevTools:Unknown Error','Unknown error - this really should not have happend')
end


%%%%%%%%%%%%%%%%%%%%Compute Controller%%%%%%%%%%%%%%%%%%%%%%%%%%


Bg=[zeros(ng,lr) B1gw B2g];
Br=[[eye(lr)*Dr;Br] zeros(lr+nwr,m+lw) ];

Lb=[Lh zeros(l+m,lr*(N-1)+nwr)]+[BX(:,1:ng)*Ag BX(:,1:ng)*B1gr BX(:,ng+1:ng+(N-1)*lr) BX(:,ng+(N-1)*lr+1:ng+N*lr)*Cr+BX(:,end-nwr+1:end)*Ar ];
Lb2=Lb(l+1:l+m,:);
Lb1=Lb(1:l,:);


%[A,B,R,L,Q]=GetSSb(SetN(P,N),gam); % for debugging only
%Xcheck=Xinfd(SetN(P,N),gam); % for debugging only
%Lbcheck=L'+B'*Xcheck*A; % for debugging only
%Lbcheck-Lb % for debugging only

Jh=[eye(l) zeros(l,m); zeros(m,l) -gam^2*eye(m)];
[vr3,dr3]=svd(Rb3);
V12=sqrt(dr3)*vr3;
if norm(V12'*V12-Rb3)>1e-6
	error('V12 Decomposition has failed')
end

'Lnab';
Lnab=Lb1-Rb2'*Rb3inv*Lb2;

'nabinv';
nabinv=inv(nab);

if ~bFI
    if P.N~=N 
        error('OF controller computation currently fails if P.N is incorrectly set')
    end
    'A';
    A;
    'B1';
    B1;
    'At';
    At=A-B1*nabinv*Lnab;
    Ct1=V12*Rb3inv*(Lb2-Rb2*nabinv*Lnab);
    Ct2=C2-D21*nabinv*Lnab;
    'Ct';
    Ct=[Ct1;Ct2];

    Bttil=[B1*inv(V21) zeros(n,m)];
    V21inv=inv(V21);
    V12inv=inv(V12);
    Dttil=[V12*Rb3inv*Rb2*V21inv  eye(m);D21*V21inv zeros(q,m)];

	St=Dttil*Jh*Dttil'+Ct(:,1:ng)*Zg*Ct(:,1:ng)';
	Mt=Bttil*Jh*Dttil'+At(:,1:ng)*Zg*Ct(:,1:ng)';

	St1=St(1:m,1:m);
	St2=St(1:m,1+m:end);
	St3=St(1+m:m+q,1+m:m+q);
	Mt2=Mt(:,m+1:m+q);
	St3inv=inv(St3);

	AK=At-Mt2*St3inv*Ct2-B2*V12inv*(Ct1-St2*St3inv*Ct2);
	BK=Mt2*St3inv-B2*V12inv*St2*St3inv;
	CK=-V12inv*Ct1+V12inv*St2*St3inv*Ct2;
	DK=-V12inv*St2*St3inv;

	CKt=inv(eye(m)+DK*D22)*CK;
	DKt=inv(eye(m)+DK*D22)*DK;
	AKt=AK-BK*D22*CKt;
	BKt=BK-BK*D22*DKt;
	K=ss(AKt,BKt,CKt,DKt,P.Ts);
else
	K=-inv(Rb3)*[Lb2 Rb2];
end

E=-gam^2*log(abs(det(-gam^-2*nab)));


%*************************************************************

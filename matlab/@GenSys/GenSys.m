function sys=GenSys(G,q,m)
% sys=GenSys(G,q,m)
% Constructs a generalised plant (i.e. one with w and u as inputs
% and z and y as outputs) from the ss system, G
% q is the number of measurements (i.e. length(y))
% m is the number of controls (i.e. length(u))

if nargin==0
    sys.n=0; 
    sys.p=0;
    sys.q=0;
    sys.l=0;
    sys.m=0;
    sys=class(sys,'GenSys',ss);
else

    [A,B,C,D]=ssdata(G);
    n=size(A,1);
    l=size(B,2)-m;
    p=size(C,1)-q;

    sys.n=n; % number of states
    sys.p=p; % number of controlled variables (i.e. length (z))
    sys.q=q; % number of measurements
    sys.l=l; % number of disturbances (i.e. length(w))
    sys.m=m; % number of controls
    sys=class(sys,'GenSys',ss(G),PublicProperties());
end

superiorto('ss')
function [Kopt,gamopt,Eopt]=MkOptHinfdK(P,gam0,h0,hmin,bFI,option)
% [Kopt,gamopt,Eopt]=MkOptHinfdK(P,N,gam0,h0,hmin,bFI)
%
% Function to compute hinf-optimal controller by iterating over gamma.
%
% gam0: initial guess at minimum gamma
% h0:   initial interval size for binary search
% hmin: min interval size for binary search
% bFI:  boolean, true if only the full information controller is required

nTriesMax=100;

bVerbose=false;
if nargin>5
    if strcmp(option,'verbose')
        bVerbose=true;
    end
end
gam=gam0;
h=h0;
%find feasible start point
nCounter=0;

while true
    nCounter=nCounter+1;
    try
        if bVerbose; disp(sprintf('Trying gam=%4.4d',gam));end
        [Kopt,Eopt]=MkHinfdK(P,gam,bFI);
		gamopt=gam;
        if bVerbose;disp('...sucess!');end
        bFail=false;
    catch
        if bVerbose;disp(['failed:' lasterr ]);end
        bFail=true;
    end
    
    if bFail==true
        gam=gam+h;    
    else
        break;
    end
    
    if nCounter>nTriesMax
        error('PrevTools:MaxTriesExceeded','Could not find gamma for which K exists')
    end
end
%now decrease gam until not feasible
while true
    if bFail==true
        h=h/2;
        if h<hmin
            break
        end
        gam=gam+h;
    else
        gam=gam-h;
    end
    
    try
        if bVerbose;disp(sprintf('Trying gam=%4.4d',gam));end
        [Kopt,Eopt]=MkHinfdK(P,gam,bFI);
        gamopt=gam;
        bFail=false;
        if bVerbose;disp('...sucess!');end
    catch
        bFail=true;
        if bVerbose;disp(['failed:' lasterr ]);end;
    end
    
end

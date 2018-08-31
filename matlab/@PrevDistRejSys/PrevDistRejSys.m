classdef PrevDistRejSys < DistRejGSys

    properties
        Wr
        Wz
        Ww
        nwr
       
        N
        G
        ng
        GW
    end
    
    methods
        
    
        function psys=PrevDistRejSys(G,N,Wz,Wr,Ww)
        % P=PrevDistRejSys(G,N,Wz,Wr,Ww)
        %
        % Constructor for PrevSys
        % N : Preview length
        % Wz : Weighting function for z (can be an lti sys or simple matrix gain)
        % Wr : Weighting function for r (can be an lti sys or simple matrix gain)
        % Ww : Weighting function for w (can be an lti sys or simple matrix gain)
        % G : Object of class DistRejGSys whose output is to be regulated

        lw=G.lw;
        lr=G.lr;
        p=G.p;
        qg=G.q;
        l=G.l;
        m=G.m;
        ng=G.n;

        

        if nargin>=3;
            if size(ss(Wz))==p
                Wz=ss(Wz);
            else
                error('Incorrect dimensions for Wz')
            end
        else
            Wz=ss(eye(p));    

        end

        if nargin>=4;
            if size(ss(Wr))==lr
                Wr=ss(Wr);
                nwr=size(Wr.A,1);
            else
                error('Incorrect dimensions for Wr')
            end
        else

            nwr=0;
            Wr=ss(eye(lr));
                  
        end

        if nargin==5;
            if size(ss(Ww))==lw
                Ww=ss(Ww);
            else
                error('Incorrect dimensions for Ww')
            end
        else
            Ww=ss(eye(lw));
    
        end

        if nargin>5 | nargin<2;
            error('Incorrect number of input arguments')
        end


        % Private method used to of generate the augmented state space representation 
        % of a preview system
        [ng,p,q,l,m]=Getsz(G);
        lr=G.lr;
        lw=G.lw;
        

        % absorb Wz and Ww into G to give GW
        Wwru=ss(eye(l+m));
        Wwru(lr+1:l,1+lr:l)=Ww;
        Wyz=ss(eye(p+q));
        Wyz(1:p,1:p)=Wz;
        GW=DistRejGSys(Wyz*G*Wwru,q,m,lr);
        
        [Ag,B1g,B2g,C1g,C2g,D11g,D12,D21g,D22g,Ts]=GetSS(GW);
        [B1gr,B1gw,D11gr,D11gw,D21gr,D21gw]=GetSSrw(GW);

        ng=size(Ag,1); % redefine ng

        nps=N*lr;

        Ap=[zeros(nps-lr,lr) eye(nps-lr,nps-lr);zeros(lr,nps)];
        Bp=[zeros(nps-lr,lr);eye(lr)];
        Cp=[eye(lr) zeros(lr,nps-lr)];

        C1=[C1g D11gr*Cp ];

        C2=[C2g D21gr*Cp;zeros(lr,ng+nps)];

        D21=[ zeros(q,lr) D21gw;  eye(lr) zeros(lr,lw) ];

        D11=[ zeros(p,lr) D11gw];

        D22=[D22g; zeros(lr,m)];

        B=[zeros(ng,lr) B1gw B2g;Bp zeros(nps,m+lw)];

        A=[Ag B1gr*Cp;zeros(nps,ng) Ap];
        C=[C1;C2];
        D=[D11 D12; D21 D22];

        Wwru=ss(eye(l+m));
        Wwru(1:lr,1:lr)=Wr;
        psys=psys@DistRejGSys(ss(A,B,C,D,G.Ts)*Wwru,q+lr,m,lr);
        psys.Wr = Wr;
        psys.Wz=Wz;
        
        psys.Ww=Ww;
        psys.nwr=nwr;
        psys.lr=lr;
        psys.lw=lw;
        psys.N=N;
        psys.G=G;
        psys.ng=ng;
        psys.GW=GW;
        
        end
    end
end

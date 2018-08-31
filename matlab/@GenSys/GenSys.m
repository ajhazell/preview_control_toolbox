classdef GenSys < ss

    properties
        n
        p
        q
        l
        m

    end
    methods
        function obj=GenSys(G,q,m)
            % obj=GenSys(G,q,m)
            % Constructs a generalised plant (i.e. one with w and u as inputs
            % and z and y as outputs) from the ss objtem, G
            % q is the number of measurements (i.e. length(y))
            % m is the number of controls (i.e. length(u))
            [A,B,C,D,Ts] = ssdata(G);
            obj = obj@ss(A,B,C,D,Ts); % no idea why we need to unpack and pack, but something breaks if we construct from ss(G)
 
            if nargin==0
                obj.n=0; 
                obj.p=0;
                obj.q=0;
                obj.l=0;
                obj.m=0;
                
            else

                [A,B,C,D]=ssdata(G);
                n=size(A,1);
                l=size(B,2)-m;
                p=size(C,1)-q;

                obj.n=n; % number of states
                obj.p=p; % number of controlled variables (i.e. length (z))
                obj.q=q; % number of measurements
                obj.l=l; % number of disturbances (i.e. length(w))
                obj.m=m; % number of controls
                
            end
            
        end
    end
end

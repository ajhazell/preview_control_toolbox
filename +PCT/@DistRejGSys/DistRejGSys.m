classdef DistRejGSys < PCT.GenSys

    properties
        lr
        lw
    end
    
    
    methods
       

        function gsys=DistRejGSys(G,q,m,lr)
            % DistRejGSys(G,q,m,lr)
        %
        % Constructor for DistRejGSys - Helper class for PrevDistRejSys
        %
        % q: : Dimension of the measurement (i.e. length(y))
        % m  : Dimension of the control (i.e. length(u))
        % lr : Dimension of reference to be tracked
        % G  : Object of class ss 

        
        gsys=gsys@PCT.GenSys(G,q,m)
        

        l=size(G.B,2)-m;
        lw=l-lr;

        if lw<0
            error('PrevTools:lwnegative','Inferred value of lw is negative')
        end

        gsys.lr=lr;
        gsys.lw=lw;

        
        end
    end
end

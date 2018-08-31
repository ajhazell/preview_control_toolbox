classdef LQRTrackSys < PCT.PrevTrackSys
% P=LQRTrackSys(G,R,Q,N)
%   
% Computes generalised plant for LQR preview tracking formulation
%
% Note that the cost function has the form J=sum [e'Qe+u'Ru]
% i.e. Q is an error weighting, not a state weighting
%
% N is the preview length

    properties
        Q
        R
    end
    
    methods
 
        function psys=LQRTrackSys(G,R,Q,N)
        % P=LQRTrackSys(G,R,Q,N)
        %   
        % Computes generalised plant for LQR preview tracking formulation
        %
        % Note that the cost function has the form J=sum [e'Qe+u'Ru]
        % i.e. Q is an error weighting, not a state weighting
        %
        % N is the preview length
        %
        
        m=size(G,2);
        p1=size(G,1);
        lr=p1;
        Gtrack=[G;eye(m)];
        We=chol(Q);
        Wu=chol(R);
        
        Wz=[We zeros(p1,m);zeros(m,p1) Wu];
        psys = psys@PCT.PrevTrackSys(PCT.GenSys(Gtrack,0,m),N,lr,Wz);

        psys.Q=Q;
        psys.R=R;

        end
    end
end


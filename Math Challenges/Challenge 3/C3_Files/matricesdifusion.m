function [M,K]=matricesdifusion(m,k,nui,nud,J)
% m k son o bien vectores de longitud J y J-1, y no se declara J,
% o bien constantes y se declrara J, en cuyo caso la funci\'on
% construyes lo vectores
% nui y nud son sendos coeficietnes de transmisi\'on en los extremos.
% nu = inf se acepta, que corresponde a la condici\'on Dirichlet.
% nu = 0 es la condici\'on Neumann.
% otros valores corresponden a considiones Robien
%
    if length(m)==1
        m=m*ones(J,1);
    end
    if length(k)==1
        k=k*ones(J-1,1);
    end
M=diag(m);
K=diag(k,-1)+diag(k,1)-diag([0;k]+[k;0]);
    if (nui==inf)
        K(1,1)=-1;
        K(1,2)=0;
    else
    K(1,1)=K(1,1)-nui;
    end
    if (nud==inf)
        K(end,end)=-1;
        K(end,end-1)=0;
    else
    K(end,end)=K(end,end)-nud;
    end
end

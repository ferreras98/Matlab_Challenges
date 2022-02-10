function [Qk,sk,Pk,Ak] = codificaxsvd(A,criterio,tolrel)
% criterio = p, con p=0, es la norma de tipo operador
%               con 1 \le p \le inf, es la norma p de X(:)
% criterio ='fro', es lo mismp que p=2
% criterop = 'op', es lo mismo que p=0;
if nnz(A)==0;
    Ak=zeros(size(A));
    Qk=[];
    sk=[];
    Pk=[];
    return
end
[Q,s,P]=mysvd(A);
if strcmp(criterio,'fro')
    cs=cumsum(s.^2);
    k=find(cs >= cs(end)*(1-tolrel)^2,1);
elseif strcmp(criterio,'op')
    k=find(s<=tolrel*s(1),1);
end
Qk=Q(:,1:k);
sk=s(1:k);
Pk=P(:,1:k);
Ak=Qk*diag(sk)*Pk'; % es opcional, si no se pide no se calcula

end


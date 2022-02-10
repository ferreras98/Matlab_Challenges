function [Q,S,P] = mysvd( A )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A=Q*diag(S)*P'
%P,Q tienen columnas ortonormales
%long columnas de Q = long columas de A
%long columnas de P = long filas de A
%S es un vector columna de elementos estrictamente positivos y
%ordenados de forma decreciente. Los elementos de sigma son los 
%valores singulares de A
% Las columnas de P son los llamados autovectores derechos de A
% Las columnas de Q son los llamados autovectores izquierdos de A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Si A es m x n y long(S)=k,
%       Q es m x k
%       P es n x k
%       A = sum_{j=1}^k S(j) * (col j de Q) * (col j de P)'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%              A x = sum_{j=1}^k S(j) * < x,p_j> q_j
%
%donde q_j y p_j son las columnas de Q y de P. Otra forma es
%
%               A = sum_{j=1}^k S(j) *  q_j*p_j'
%
%que descompone A en una suma de matrices de rango 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nnz(A)==0 || isempty(A))
    
    S=[];
    P=[];
    Q=[];
    
else
[m,n]=size(A);

    if (m>=n)
      B=A'*A;
      [P,d]=eig(-B,'vector');
      d=-d;
      r=length(d(d>0));
      S=d(1:r).^0.5;
      P=P(:,1:r);
      Q=(A*P).*repmat(1./S',m,1);
%       Q=A*P;
%          for k=1:r
%          Q(:,k)=(1/S(k))*Q(:,k);
%          end
      
    else
       B=A*A';
       [Q,d]=eig(-B,'vector');
       d=-d;
       r=length(d(d>0));
       S=d(1:r).^0.5;
       Q=Q(:,1:r);
       P=(A'*Q).*repmat(1./S',n,1);
%        P=(A')*Q;
%           for k=1:r
%           P(:,k)=(1/S(k))*P(:,k);
%           end
     end
    
end 


end


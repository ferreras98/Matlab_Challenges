function CL = legendre(N)
% CL es una matriz (N+1) \times (N+1).
% La columna CL(:,k), 1\le k \le N+1, da los coeficientes del polinomio de 
% Legendre L_k(x) (normalizaci\'on L_k(1)=1) para el intervalo [-1,1].
% Los coeficientes siguen el orden creciente de las potencias:
%
% L_k(x) = \sum_{j=1}^{k+1} CL(j,k+1) x.^{j-1}, \qquad 1 \le j \le k+1.
% 
% La funci\'on se basa en la recurrencia de Bonnet:
%
% L_0=1;
% L_1=x;
% (2n+1)xL_n(x)=(n+1)L_{n+1}(x)+nL_{n-1}(x)
%
CL=zeros(N+1,N+1);
CL(1,1)=1;
CL(2,2)=1;
for n=1:N-1
    pn=CL(:,n+1);
    pnm1=CL(:,n);
    CL(:,n+2)=((2*n+1)*circshift(pn,1)-n*pnm1)/(n+1);
end

end


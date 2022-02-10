function [x,defect,niter] = newton(f,Df,x0,nmax,TOL,hinJ )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if ~isa(Df,'function_handle')
    J=@(x,y) jacobianaaprox(f,x,y,hinJ);
else
    J=@(x,y) Df(x);
end
xold=x0;
niter=0;
yold=f(xold);

xnew=xold-J(xold,yold)\yold;
niter=niter+1;
defect=norm(xnew-xold);

while  ((niter < nmax) && (defect > TOL))
    
    xold=xnew;
    yold=f(xold);
    xnew=xold-J(xold,yold)\yold;
    niter=niter+1;
    defect=norm(xnew-xold);
    
end

x=xnew;

end


function J=jacobianaaprox(f,x,y,h)
N=length(x);
I=eye(N);
J=zeros(N);
for n=1:N
    en=I(:,n);
    J(:,n)=f(x+h*en)-y;
end
J=(1/h)*J;
end

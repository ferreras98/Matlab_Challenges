function nodos = chebnodes(J,a,b)
% (1) Para a=-1, b=1: Se trata de los valores cos(\theta) tales que
%     cos(J\theta)=0.
%
%         cos(J \theta)=0 \Rigtharrow J\theta=\pi/2+m*\pi=
%              \Rigtharrow \theta= (\pi/2+m*pi)/J
%
% Tomamos los ceros en (0,pi), sin repeticiones
%
%              \theta_m=(\pi/J)*(1/2+m), \qquad m=0:(J-1)
%
% lo que lleva al vector de nodos (de  Chebyshev) 
%
%              nodos=cos((pi/J)*(1/2+(0:J-1))';
%
% Notemos que
%              nodos(m) = real(exp(1j*pi/N)*z_J^m), 
%
% siendo N=2*J y z_J=exp(2*pi*m/J),  la ra\'{\i} primitiva J-\'esima de la 
% unidad.

% (2) Para (a,b), llamamos R=(b-a)/2, c=(a+b)/2,
%
%              nodos=c+R*cos((pi/J)*(1/2+(0:J-1))';
%
%  Generamos el vector nodos en columna; length(nodos)=J y est\'an 
%  ordenados de forma decrecente.

% Por defecto, declarando s\'olo chebnodes(J), adopamos a=-1,b=1;

nodos=cos((pi/J)*(1/2+(0:J-1)))';

    if (nargin==3)
        R=(b-a)/2;
        c=(a+b)/2;
        nodos=c+R*nodos;
    end

end

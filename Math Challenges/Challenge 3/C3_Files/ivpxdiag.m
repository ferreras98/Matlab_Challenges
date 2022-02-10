function U = ivpxdiag(varargin)
%%%%%%%%%%%%%%%% FORMATO Y PROP\'OSITO DE LA FUNCI\'ON %%%%%%%%%%
% Formatos de entrada:
% (K)
% (K,true/false), 
% (M,K)
% (M,K,true/false)
%
% Se supone que M es regular y que la matriz A=M^{-1}K es diagonalizable. 
% La funci\'on ivpxdiag resuelve el IVP.
%
% (IVP)   Mu'(t) = K u(t), \qquad u(0)= u0.
%
% Notemos que el tiempo inicial en (IVP) se toma 0. Para otro tiempo 
% t0 basta con trabajar con t-t0.
%
% Si no se declara M, entonces M=Id. Si M es diagonal, podemos emplear
% el vector diag(M).
% 
% El argumento l\'ogico se declara true cuando
%     (1) M es definida positiva y adem\'as
%     (2) K es herm\'{\i}tica o antiherm\'{\i}tica.
% Por defecto es falso, cosa que tambi\'en se puede declarar o no. El caso
% true conlleva una simplificaci\'on y una mejora de la precisi\'on, que
% conviene aprovechar. 

% La funci\'on ivpxtrig va a resolver el (IVP) generando el 
% function_handle U de la soluci\'on como funci\'on de t y de u0. 
% El output U es pues una funci\'on de dos variables:
%     
% U(t,u0)= valor en t de la soluci\'on de (IVP) con dato inicial u0 en
%           el tiempo inicial t_0=0.
%
% Cuando t sea un vector de tiempos y u0 sea fijo, U(t,u0) es una matriz: 
% la columna n de U(t,u0) es el valor de la soluci\'on de (IVP) en tiempo t(n).
%
%%%%%%%%%%%%%%%%%% FUNDAMENTO MATEM\'ATICO %%%%%%%%%%%%%%%%%%%%%%%
%
% A = M^{-1}K tiene cierta dimensi\'on N \times N. Formalmente, la 
% soluci\'on de (IVP) es
%
%      u(t) = \expm(tA) u0,
%
% representaci\'on que afirma dos cosas muy importantes: (1) que la
% soluci\'on es \'unica y (2) que depende linealmente del dato inicial u0.
%
% Si u0=p es un atovector de A, con autovalor d, entonces
%
%   \expm(tA) p = \exp(td)p   (modo propio: la soluci\'on es una
%   modulaci\'on en tiempo de p, con el factor exponencial indicado).
%
% [P,D]=eig(A) nos da una base autovectores (las N columnas de P) y los
% autovalores correspondientes (la diagonal de la matriz diagonal D).
% Llamemos p_n =P(:,n),  1 \le n \le N, y d=diag(D), de forma que
%
%  A*P = P*D , o bien, Ap_n = d(n) p_n, 1 \le n \le N.
%
% Sabiendo ya que los N modos propios
%
%      \expm(tA) p_n = \exp(t d(n)) p_n 
%
% forman una base del espacio de soluciones de (IVP), dado u0 necesiamos
% encontrar un vector c0 tal que
%
% (An\'alisis)     u0 = \sum_{n=1}^N c0(n) p_n, (basta tomar c0=P\c0).
%
% Encontrado c0, la soluci\'on de (IVP) es
%
%   u(t)=\sum_{n=1}^N c0(n)\exp(td(n)) p_n, o bien,
%
% (S\'{\i}ntesis)         u(t)=P*diag(\exp(t*d).*c0)  (empleamos vectores,
%                                                      y luego ponemos en
%                                                      la diagonal).
%                                                       
%
% Una situaci\'on muy especial e importante se presente cuando M es definida 
% positiva y  K es herm\'{\i}tica o antiherm\'{\i}tica.
% La diagonalizaci\'on de A puede lograrse entones de una forma muy estable 
% mediante
%
%                        [Q,E]=eig(K,M)
%
% Se plantea el problema de autovectores del problema generalizado 
%
%            K v = \lambda M v. 
%
% Los autovalores se guardan en e=diag(E) y los autovectores en las
% columnas de Q. De este modo se va a cumplir
%
%        K*Q(:,n) = e(n)*M*Q(:,n), o bien K*Q=M*Q*E.
%
% En las condiciones se\~naladas siempre posible realiar esta
% diagonalizaci\'on, seleccionando adem\'as los autovectores para que
% formen un base ortonormal respecto del producto interno definido por M:
%
%             transpose(Q)*M*conj(Q) = Id.
% 
% En el lenguaje de A=M^{-1}K tendremos
%
%      A*Q=Q*E
% 
% de forma que Q es una base de autovectores para A y E es diagonal con los
% autovalores. Si ejecutamos [P,D]=eig(A), podemos encontramos otra base de 
% autovectores y los autvalores en otro orden. Entendido esto, trabajamos
% directamente con 
%
%     [P,D]=eig(K,M).
%
% y significado de P y de D es el mismo (diagonalizaci\'on de A=M^{-1}K).
% Como ahora
%
%      P^{-1}=P'*conj(M)
%
% resulta muy sencillo (y num\'ericamente muy preciso) resolver P*c0=u0
% mediante
%
% (An\'alisis)       c0=P'*(conj(M)*u0).
%    
% Adem\'as, muy importante, esta idea admite una extensi\'on natural para 
% tratar problemas en los que ya no es adecuado pensar en el formato matriz 
% por vector.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nunca sobra mirar si al menos el n\'umero de argumentos es el esperado.
narginchk(1,3);
if (nargout~=0)
nargoutchk(1,1);
end
% Comenzamos leyendo varargin y clasificando (caso 1: M = Id; caso 2: M
% es diagonal; caso 3: M is full; status=true lee/completa el l\'ogico )
%
if ( (nargin==1) || (nargin==2 && islogical(varargin{2}))  )
    K=varargin{1};
    status=false;
    caso=1;
    if (nargin==2)
        status=varargin{2};
    end
else
    M=varargin{1};
    K=varargin{2};
    [m,n]=size(M);
    caso=3;
        if (m~=n)
            caso=2;
        end
    status=false;
        if (nargin==3)
            status=varargin{3};
        end
end
% Al final nos vamos a preguntar si hay que tomar parte real. Por el
% momento vamos tomando nota
if (caso==1)               
   outinreal=isreal(K);    
else
   outinreal = (isreal(M) && isreal(K));
end
% Pasamos a formar y diagonalizar A y, fundamental, a definir analysis
if ~status                         
    if ( caso==1 )
        [P,D]=eig(K);
    elseif (caso==2)
        [P,D]=eig(diag(1./M)*K);
    else
        [P,D]=eig((M\eye(size(M)))*K);
    end
    analysis=@(u0) P\u0; %  Ser\'{\i}a ideal el poder emplear
                         %       [LP,UP,p,invq]=mylu(P);
                         %       analyze=@(u0) lusolver(LP,UP,p,invq,u0); 
                         %  si estas funciones estuvieran implementadas en 
                         %  MATLAB nativo, algo no disponible desde hace
                         %  varias versiones (nuestra implementaci\'on de
                         %  tales funciones es impecable matem\'aticamente
                         %  hablando, pero mylu es muy lenta a patir de 
                         %  length=100).
                         % 
                         %  Una alternativa, s\'olo cuando haya   
                         %  que resolver despu\'es muchos sistemas (con los
                         %  mismos coeficientes, pero con distintos datos 
                         %  iniciales) puede ser
                         %         inversaP=P\eye(size(P))
                         %         analyze=@(u0) inversaP*u0;
else
    if (caso==1)
        [P,D]=eig(K);
        analysis=@(u0) P'\u0;
    elseif (caso==2)     
        [P,D]=eig(K,diag(M));
        analysis=@(u0) P'*(conj(M(:)).*u0);
    else
        [P,D]=eig(K,M);
        analysis=@(u0) P'*(conj(M)*u0);
    end
end
% Ahora podr\'{\i}amos clear K M.
d=diag(D);
% Para retroceder a las coordenadas de origen
synthesis=@(C) P*C;
% Ahora podr\'{\i}amos clear D y P.

% Presentamos directamente la soluci\'on en forma de function_handle.
% El trabajo se realiza en solucion y ahora simplemente recogemos:
U=@(t,u0) solucion(d,t,u0,analysis,synthesis,outinreal);
end


function S=solucion(d,t,u0,analysis,synthesis,outinreal)
% Este programa es un ejemplo fundamental para entender multitud
% de situaciones del mundo lineal.

% Analizamos el dato:
c0=analysis(u0(:)); %Hay que forzar el formato columna en u0.

% Resolvemos en las componentes transformadas:
S=zeros(length(c0),length(t));
    for n=1:length(t)
        S(:,n)=exp(t(n)*d).*c0; % Simplificamos mucho con .*  !!
    end
    
% Sintetizamos el resultado, sobreescrbiendo en S:
S=synthesis(S);

% Por \'ultimo, nos preguntamos si tomamos la parte real
    if (outinreal && isreal(t) && isreal(u0))
        S=real(S);
    end
end




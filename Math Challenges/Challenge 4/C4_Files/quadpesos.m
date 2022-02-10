function [w,m,defect] = quadpesos(theta,tol)
% theta es un vector de N+1 abscisas sin dimensiones, los nodos
% de la cuadratura vistos en [0,1]. 
% w es un vector del mismmo tipo que theta que  nos da los pesos
% de la cuadratura con nodos theta y con grado de precisi\'on 
% al menos s-1.
% m da el grado de exactitud de la f\'ormula hasta tol.
% Por defecto tol =10^(-16).
% defect es el error para el grado m+1.
% Si defect es muy peque\~no indica que quiz\'as el orden fuera 
% al menos m+1; en todo caso, a efectos pr\'acticos la f\'ormula
% se comporta como una grado de precisi\'on \le m+1.
%%%%
% Para la matriz de Vandermonde hay que tratar theta como fila.
% La soluci\'on sale en columna y al final damos el formato que 
% tuviera theta en la entrada.
%%%

if (nargin==1)
    tol=10^(-16);
end
[p,q]=size(theta); % el la siguiente l\'{\i}nea se puede perder
                   % el formato de theta.
theta=theta(:)'; % lo vamos a pensar en fila.
N=length(theta)-1;% es el N de la teor\'{\i}a.
F=repmat(theta,N+1,1);% o  bien ones(N+1,1)*theta
P=repmat((0:N)',1,N+1); % o bien (0:N)'*ones(1,N+1);
V=F.^P; % matriz de Vandermonde modelada con la fila theta
Iexact=1./(1:N+1)'; % integrales de los monomios en [0,1],
                    % desde la potencia 0 a la N.
w=V\Iexact;% pesos en columna para alcanzar al menos grado de p
           % precisi\'on N.
                            
% De momento el grado m es al menos N. Probamos si es mayor:

m=N;
trym=m+1;
exact=1/(trym+1); %vemos si es exacta para la porencia trym
quad=(theta.^trym)*w;% recordemos que por ahora theta es fila
                     % y w es columna.
defect=abs(exact-quad);
while (defect < tol)
    m=m+1;
    trym=m+1;
    exact=1/(trym+1);
    quad=(theta.^trym)*w;
    defect=abs(exact-quad);
end
w=reshape(w(:),p,q); % para dar a w el formato que tiene theta
                     % a la entrada.

end


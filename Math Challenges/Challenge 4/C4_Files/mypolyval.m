function y = mypolyval(varargin)
%
% varargin=[coeficientes, centros, x}
%
%Se trata de evaluar un polinomio p_N(x) de grado a lo sumo N.
%La entrada de la funcion puede ser:
%%%%%
%a) Polinomio en potencias de x: Si p_N(x) = c_0 + c_1 x + ...+c_N x^n
% en cuyo caso varargin == {coef,x}, siendo coef el vector que guarda los
% coeficientes, coef(k)=c_{k-1}, 1 \le k \le N+1. A diferencia de MATLAB,
% los coeficientes se ordenan en potencias CRECIENTES de x.
%
% b) Polinomio en potencias de (x-x_0): Si 
% p_N(x) = c_0 + c_1 (x-x_0) + ...+c_N (x-x_0)^n, haremos
% varargin=  {coef,x_0, x}, donce coef es como en a)
%
% c) Polinomio en base de Newton:
% p_N(x) = c_0 + c_1(x-x_0)+ c_2(x-x_0)(x-x_1) + ... 
%               + c_N(x-x_0)(x-x_1)...(x-x_{N-1}),
% haremos {coef,centros, x}, donde centros es el vector [x_0, ...,x_{N-1}].
%
% Tanto coef como centros se pueden manejar en formato fila o columna.
% Adem'as acepta un vector de centros que contenga componentes 
% extra.
% La salida y tiene el formato de x, que pueden ser vectores e 
% incluso matrices. La componente j de y es la evaluacion de p_N en la componente
% j de x.
%
% Se utiliza el algoritmo de Horner (como la regla de Ruffini).
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargin==3)
    coef=varargin{1};
    centros=varargin{2};
    x=varargin{3};
    if (length(centros)==1)
        centros=centros*ones(size(x));
    else
        centros=centros(1:(length(coef)-1));
    end
end

if (nargin==2)
    coef=varargin{1};
    x=varargin{2};
    centros=zeros( (length(coef)-1),1);
end

Nm1=length(coef);
y=coef(Nm1);%comienza siendo numero
for n=(Nm1-1):(-1):1
    y=coef(n)+y.*(x-centros(n));%tras el primer paso y se convierte 
                                %en array del mismo tipo que x
end
    

end


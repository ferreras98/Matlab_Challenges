%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 1 - M�TODOS NUM�RICOS                                             %
%% Parte 1 - Una funci�n MATLAB                                           %
%%                                                                        %
%% Antonio Ferreras Extremo                                               %
%% Diciembre 2019                                                         %
%%                                                                        %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = amprespfrec(p, qi, qd, qs, w, F)

% Construimos la matriz tridiagonal
d = 1i*w*p + qd;    % Diagonal Principal
a = qi;             % Subdiagonal inferior
b = qs;             % Subdiagonal superior

% Matriz LU or el m�todo de Thomas
[li,ud,us] = luthomas(a,d,b);

% Inversa y c�lculo de la soluci�n
U = solverthom(li,ud,us,F);

end
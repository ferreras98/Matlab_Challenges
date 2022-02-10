%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 1 - MÉTODOS NUMÉRICOS                                             %
%% Parte 1 - Una función MATLAB                                           %
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

% Matriz LU or el método de Thomas
[li,ud,us] = luthomas(a,d,b);

% Inversa y cálculo de la solución
U = solverthom(li,ud,us,F);

end
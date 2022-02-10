%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 3 - Asignatura: MÉTODOS NUMÉRICOS                                 %
%% Parte 2: Construcción de la función E = ivpleftrieght(ML, KL, MR, KR)  %
%% Antonio Ferreras Extremo                                               %
%% Diciembre 2019                                                         %
%%                                                                        %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% La salida es la función handle, E = E(ts, W0) que resuelve:             %
%      W'(t) = L*W(t) + W(t)*R sujeta a W(0) = W0                         %
% para los tiempos definidos en ts, siendo:                               %
%    -  L = ML^{-1} * KL                                                  %
%    -  R* = MR^{-1} * KR                                                 %
% Entradas:                                                               %
%       Matrices: ML, KL, MR, KR                                          % 
% Salida: Función handle E(ts, W0, tol)                                   %
%           - ts: Vector de tiempos                                       %
%           - W0: Matriz mxn de condiciones inciales                      %
%           - tol: tolerancia relativa de la solución                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = ivpleftright(ML, KL, MR, KR)

    V = ivpxdiag(ML, KL);
    U = ivpxdiag(MR, KR);

    % Devolvemos el handler de la función
    E = @(ts, W0, tol) solucion(ts, W0, tol, V, U);
end

function W=solucion(ts, W0, tol, V, U)
    
    % Distinguimos si la tolerancia es idénticamente 0 (al elevar al cuadrado 
    % se pierden términos.
    if (tol == 0)
        [Q, S, P] = mysvd(W0);
    else
        [Q, S, P] = codificaxsvd(W0, 'fro', tol);
    end

    % Número de autovalores retenidos
    r = length(S);
    % Instantes de tiempo dónde calculamos la solución
    n_t = length(ts);

    % Reservamos el espacio de la solución
    W = zeros([size(W0), n_t]);

    for k=1:r
        v0 = Q(:,k);
        u0 = P(:,k);
        s0 = S(k);
        vt = V(ts, v0);
        ut = U(ts, u0);
        for t=1:n_t
            W(:,:,t) = W(:,:,t) + s0*vt(:,t)*ut(:,t)';
        end
    end    

end


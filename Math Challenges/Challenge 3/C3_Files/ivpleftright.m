%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 3 - Asignatura: M�TODOS NUM�RICOS                                 %
%% Parte 2: Construcci�n de la funci�n E = ivpleftrieght(ML, KL, MR, KR)  %
%% Antonio Ferreras Extremo                                               %
%% Diciembre 2019                                                         %
%%                                                                        %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% La salida es la funci�n handle, E = E(ts, W0) que resuelve:             %
%      W'(t) = L*W(t) + W(t)*R sujeta a W(0) = W0                         %
% para los tiempos definidos en ts, siendo:                               %
%    -  L = ML^{-1} * KL                                                  %
%    -  R* = MR^{-1} * KR                                                 %
% Entradas:                                                               %
%       Matrices: ML, KL, MR, KR                                          % 
% Salida: Funci�n handle E(ts, W0, tol)                                   %
%           - ts: Vector de tiempos                                       %
%           - W0: Matriz mxn de condiciones inciales                      %
%           - tol: tolerancia relativa de la soluci�n                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E = ivpleftright(ML, KL, MR, KR)

    V = ivpxdiag(ML, KL);
    U = ivpxdiag(MR, KR);

    % Devolvemos el handler de la funci�n
    E = @(ts, W0, tol) solucion(ts, W0, tol, V, U);
end

function W=solucion(ts, W0, tol, V, U)
    
    % Distinguimos si la tolerancia es id�nticamente 0 (al elevar al cuadrado 
    % se pierden t�rminos.
    if (tol == 0)
        [Q, S, P] = mysvd(W0);
    else
        [Q, S, P] = codificaxsvd(W0, 'fro', tol);
    end

    % N�mero de autovalores retenidos
    r = length(S);
    % Instantes de tiempo d�nde calculamos la soluci�n
    n_t = length(ts);

    % Reservamos el espacio de la soluci�n
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


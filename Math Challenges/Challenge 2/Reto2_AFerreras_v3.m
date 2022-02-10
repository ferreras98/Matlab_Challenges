%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 2 - Asignatura: MÉTODOS NUMÉRICOS                                             %
%%                                                                        %
%% Antonio Ferreras Extremo                                               %
%% Noviembre 2019                                                         %
%%                                                                        %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Datos del problema de la red LRC
nc = 100;        % Número de nodos
c  = 0.0001;     % Valor de la capacidad
l  = 0.01;       % Valor de la inductancia
r  = 0.05;       % Valor de la resistencia

% Cargamos los datos de voltaje medido en el último nodo
load("datosreto2.mat")

% Construimos nuestra matriz LRC (siguiendo el módulo matrizLRC.m)
C = c*ones(nc, 1);
L = l*ones(nc-1, 1);
R = r*ones(nc-1, 1);
E = -eye(nc) + diag(ones(nc-1,1), -1);
E = E(:, 1:nc-1);
A = [-diag(R./L), -diag(1./L)*E'; diag(1./C)*E, zeros(nc)];

% Cálculo de los autovectores y autovalores de la matriz
[P, D] = eig(A);

% Calculamos el valor de las autofunciones para los ts del fichero
d = diag(D).';         % La transponemos sin que calcule la conjugada!
VALORES = exp(ts(:)*d) .* P(end,:);  % Matriz de condiciones de voltaje

% Añadimos las condiciones adicionales de contorno adicionales de I(t=0) = 0
VALORES = [P(1:nc-1,:); VALORES];

% Vector de coeficientes
b = [zeros(nc-1,1); v(:)];

% Resolvemos el ajuste de mínimos cuadrados
c0 = (VALORES'*VALORES)\(VALORES'*b); 

% Calculamos las corrientes y potenciales iniciales pedidos
u0 = real(P*c0);
voltajes = u0(nc:end);

% Salida de resultados. Consideramos también aquellos < -120 V. 
solucion = find(abs(voltajes)>120);

% Salida de resultados
str1=sprintf('%d ', solucion);
str2= sprintf('%6.2f ',voltajes(solucion));
handle = msgbox( ...
    {"Condensadores que adquirieron una tensión > 120 V"; ...
     strcat(['    NODOS :' '  ' str1]); ...  
     strcat(['VOLTAJES :' '  ' str2])}, ...  
     'Reto 2: RESPUESTA');
uiwait(handle);

% Gráfica con todas las voltajes (las corrientes iniciales son 0)
handle = figure('Name', ...
    'Métodos Numéricos en Computación - Antonio Ferreras', ...
    'Numbertitle', 'off', 'MenuBar', 'none', 'Resize', 'off');
stem((1:nc), voltajes, 'g');
hold on;
title("Reto2: Voltajes en los nodos para t=0"); 
grid on; yticks([0 50 100 120 150]);
xlabel('Número de Nodo'); ylabel('Voltaje'); 
xlim([0,nc]); ylim([-5, 170]);
stem(solucion, voltajes(solucion),'b');   
plot([0,100],[120,120], 'r');
legend("V < 120 v", "V > 120 v", "");
pause(); close(handle); clear all;

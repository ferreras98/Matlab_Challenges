%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 3 - Asignatura: MÉTODOS NUMÉRICOS                                 %
%%                                                                        %
%% Antonio Ferreras Extremo                                               %
%% Diciembre 2019                                                         %
%%                                                                        %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Datos del problema
m = 450;
n = 1000;
tolrel = 0.01;
ts = (0: 10: 1000);

load('datosreto3.mat');

% Matrices de difusión bidimensionales según el enunciado del problema
[ML, KL] = matricesdifusion(1, 200, 0, 0, m);
[MR, KR] = matricesdifusion(1, 200, 0, 0, n);

% Construimos la función que resuelve ese problema de difusión
E = ivpleftright(ML, KL, MR, KR);

% Resolución del problema para los datos iniciales de 'datosreto2.mat'
W = E(ts, W0, tolrel);

% Salida de resultados de forma gráfica
mymovie2D(W, ts, 0.24);

% Dejamos limpio el entorno
clear all;

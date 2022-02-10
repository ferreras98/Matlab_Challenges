%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 3 - Asignatura: M�TODOS NUM�RICOS                                 %
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

% Matrices de difusi�n bidimensionales seg�n el enunciado del problema
[ML, KL] = matricesdifusion(1, 200, 0, 0, m);
[MR, KR] = matricesdifusion(1, 200, 0, 0, n);

% Construimos la funci�n que resuelve ese problema de difusi�n
E = ivpleftright(ML, KL, MR, KR);

% Resoluci�n del problema para los datos iniciales de 'datosreto2.mat'
W = E(ts, W0, tolrel);

% Salida de resultados de forma gr�fica
mymovie2D(W, ts, 0.24);

% Dejamos limpio el entorno
clear all;

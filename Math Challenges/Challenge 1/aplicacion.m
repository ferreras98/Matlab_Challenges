%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 1 - MÉTODOS NUMÉRICOS                                             %
%% Parte 2 - Una ap licación                                               %
%%                                                                        %
%% Antonio Ferreras Extremo                                               %
%% Diciembre 2019                                                         %
%%                                                                        %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Datos de la red
J = 10000; c = 10^-7; k = 100;

% Condiciones de contorno
kizq = 11; Hizq = 18; Vizq = 23; wizq = 12;
kder = 20; Hder = 25; Vder = -14i; wder = 12*sqrt(3);

% Para la Matriz P
p = c*ones(1,J);

% Para la Matriz Q
qd = 2*k*ones(1,J); qd(1) = qd(1) - k + kizq; qd(J) = qd(J) - k + kder; 
qi = -k*ones(1,J-1); qs = -k*ones(1,J-1);

handle = figure('Name', ...
    'Reto 1: Métodos Numéricos en Computación - Antonio Ferreras', ...
    'Numbertitle', 'off', 'MenuBar', 'none', 'Resize', 'off');
set(gcf, 'Units', 'normalized', 'outerposition', [0 0 1 1]);
sgtitle('Módulo y Fase de la solución para las diferentes frecuencias');

% Datos en forma vectorial para utilizar un bucle de cálculo
w = [0; wizq; wder];  % En columna!
vizq = [Hizq, Vizq, 0];
vder = [Hder, 0, Vder];
U = zeros(J+2, 3, 'double');
maximo = max([abs(vizq), abs(vder)]);
G = zeros(1,J);
x = (0:J+1);

for i=(1:3)  % Para cada valor de frecuencia...
    % Excitación F (de 1 a J son 0)
    G(1) = vizq(i)*kizq; G(J) = vder(i)*kder;
    U(1,i) = vizq(i);
    % Cálculo matricial
    U(2:J+1,i) = amprespfrec(p, qi, qd, qs, w(i), G);    
    U(J+2,i) = vder(i);

    % Dibujamos el módulo del fasor
    subplot(2,3,i);
    plot(0, abs(vizq(i)), 'or', J+2, abs(vder(i)), 'or', ...
        x(2:J+1), abs(U(2:J+1,i)), 'b');
    title(strcat('w = ', 32, num2str(w(i),4))); grid on;
    xlabel('J'); ylabel('Módulo - |U|'); xlim([0,J+1]); ylim([0, maximo]);

    % Dibujamos la fase del fasor
    subplot(2,3,i+3);
    plot(0, angle(vizq(i)), 'or', J+2, angle(vder(i)), 'or', ... 
        x(2:J+1), angle(U(2:J+1,i)), 'g');
    title(strcat('w = ', 32, num2str(w(i),4))); grid on;
    set(gca,'YTick',-pi:pi/2:pi) 
    set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
    xlabel('J'); ylabel('Fase - \phi(U)'); xlim([0,J+1]); ylim([-pi, pi]);
end
pause(); close(handle);

% Representación dinámica
% NOTA: La frecuencia es más lenta (no es real) para que se pueda ver
handle = figure('Name', ...
    'Reto 1: Métodos Numéricos en Computación - Antonio Ferreras', ...
    'Numbertitle', 'off', 'MenuBar', 'none', 'Resize', 'off');

npuntos = 100;  % No es neceario representar los 10.000 puntos
x=round(linspace(1,J,npuntos));

% Cálculo del tamaño de los ejes verticales
maximo = 1.05*max([abs(Hizq) + abs(Vizq), abs(Hder) + abs(Vder)]) + 1;
minimo = 1.05*min([abs(Hizq) - abs(Vizq), abs(Hder) - abs(Vder)]) - 1;

for t=(0:0.01:10)
  yizq = real(vizq*exp(1i*w*t));
  y = real(U(x,:)*exp(1i*w*t));
  yder = real(vder*exp(1i*w*t));

  plot(0, yizq, 'or', J+1, yder, 'or', x, y, 'b');
  title(strcat('Potencial para t =', 32, num2str(t,'%.2f'), 32,...
      '(Control-C para terminar)'));
  xlabel('J'); ylabel('u(t)'); grid on; 
  xlim([0,J+1]); ylim([minimo, maximo]);
  pause(0.02);
end
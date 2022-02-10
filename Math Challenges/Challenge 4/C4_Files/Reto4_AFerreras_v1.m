%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                                                        %
%% Reto 4 - Asignatura: MÉTODOS NUMÉRICOS                                 %
%%                                                                        %
%% Antonio Ferreras Extremo                                               %
%% Diciembre 2019                                                         %
%%                                                                        %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Valor de la integral (chequeo)
% https://www.wolframalpha.com/input/?i=integrate+%281%2Bx%29*cos%28sqrt%282%2Bx%5E2%29%29+dx+from+x%3D0+to+60.367

handle = figure('Name', ...
    'Reto 4: Métodos Numéricos en Computación - Antonio Ferreras', ...
    'Numbertitle', 'off', 'MenuBar', 'none', 'Resize', 'off');

%% a. CALCULO DE LOS POLINOMIOS DE LEGRENDE, CEROS, PESOS
x = 0.:0.005:1;
Nmax = 9;   % Número de polinomios de Legendre

for N = 0 : Nmax;
    CL = legendre(N+1);                        % Cálculo de coeficientes
    funcion_L = @(x) mypolyval(CL(:,end),x);   % Handler a la funcion poli.
    r0 = chebnodes(N+1);                       % Aproximacion a los ceros
    [r,~,~] = newton(funcion_L, 'aprox', r0, 5000, 1.e-15, 1./2000);
    ceros = (1+r)/2;                    % Desplazo ceros al intervalo [0 1]
    [w,~,~] = quadpesos(ceros, 1.e-14); % Cálculo de los pesos 

    % Dibujamos cada uno de los polinomios calculados"
    sgtitle(handle, strcat('Polinomio de Legrende N+1=', num2str(N+1)));
    h1 =plot(x, funcion_L(-1+2*x), 'r');     % Dibujamos la función
    hold on;
    h2 = plot(ceros, zeros(length(r)), 'go'); % Posición de los ceros
    h3 = plot(ceros, w, 'bo');                % Peso de cada cero
    hold off;    
    xlim([-0.01 1.01]); ylim([-1.05 1.05]); grid on; box on;
    xlabel('x'); ylabel('$$ L_{N+1}(x) $$', 'interpreter','latex');
    legend([h1(1) h2(1) h3(1)], {'Polinomio', 'Ceros','Pesos'}, ...
        'location', 'southeast', 'EdgeColor', 'k', 'Box', 'on');
    pause(1); 
end

hmsg=msgbox('Pulse para continuar', 'Reto 4', 'help');
set(hmsg, 'Position', [300 200 200 65]); uiwait(hmsg);
figure(handle);

%% FUNCION A INTEGRAR. REPRESENTACION
a = 0;
b = 60.367;
x = a:0.005:b;
funcion_f = @(x) (1+x).*cos(sqrt(2+x.^2));
inty = cumtrapz(x,funcion_f(x));

sgtitle(handle, '$$ f(x) = (1+x) \cos(\sqrt{1+ x^2}) $$', 'interpreter','latex');
plot(x, inty, 'b'); 
hold on
patch([x fliplr(x)], [zeros(size(x)) fliplr(inty)], 'b')
hold off; grid on; box on; 
xlim([a-2 b+2]); ylim([-65 65]);xlabel('x'); ylabel('f(x)');
pause(1);

%% b. INTEGRACION POR CUADRATURA GAUSSIANA (GAUS - LEGRENDE)

% CALCULO DE CEROS Y PESOS. Proceso similar al anterior
N = 4;
CL = legendre(N+1);
funcion_L = @(x) mypolyval(CL(:,end),x);
r0 = chebnodes(N+1);
[r,~,~] = newton(funcion_L, 'aprox', r0, 5000, 1.e-15, 1./2000);
ceros = (1+r)/2;
[pesos,~,~] = quadpesos(ceros,1.e-14);

% b.1 Regla compuesta Q[J=1:100]
JMAX = 100;             % Número de subintervalos
Q = zeros(JMAX,1);      % Reserva de espacios

% Cálculo de la integral por cuadratura
for J=1:100
    limites = linspace(a,b,J+1);   % Dividimos el intervalo
    L = (b-a)/J;                  % Tamaño de cada subintervalo
    % Acumulamos las sumas por cada subintervalo
    for j=1:J
        a0 = limites(j);
        Q(J) = Q(J) + L*sum(pesos.*funcion_f(a0+ceros*L));
    end
end

hmsg = msgbox(strcat('Q[70] =', {' '}, num2str(Q(70), 16)),...
    'Valor Integral', 'help');
set(hmsg, 'Position', [300 200 200 65]);
uiwait(hmsg);
figure(handle);

sgtitle('Convergencia de Q(J)');
plot(1:JMAX,Q, 'b-o');
xlabel('J'); ylabel('Q[J]'); grid on; pause(1);

hmsg=msgbox('Pulse para continuar', 'Reto 4', 'help');
set(hmsg, 'Position', [300 200 200 65]);
uiwait(hmsg);

% b.2 Plot logarítmico de las diferencias
figure(handle);
sgtitle('Error en el cálculo de cuadratura');
loglog(1:69, abs(Q(1:69)-Q(70)));   % Ecala de la pendiente
hold on;
loglog([15 60], [abs(Q(15)-Q(70)) abs(Q(60)-Q(70))], 'r-*'); hold off;
xlabel("J"); ylabel('$$ |Q[J]-Q[70]| $$', 'interpreter','latex');

% Cálculo y plot de la pendiente
pendiente = (log10(abs(Q(60) - Q(70))) - log10(abs(Q(15) - Q(70)))) ... 
           /(log10(60) - log10(15));

hmsg = msgbox(strcat('Valor Pendiente :', {' '}, num2str(pendiente)),...
    'Reto 4: Pendiente 14 < J < 61', 'help');
set(hmsg, 'Position', [300 200 200 65]); uiwait(hmsg); delete(handle);
clear;
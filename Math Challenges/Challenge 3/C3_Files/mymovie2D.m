function  mymovie2D(varargin)
% varargin ={U,ts}, {U,a,b,c,d,ts}, {U,ts,espera},{U,a,b,c,d,ts,espera}
%
% ts es un vector de tiempos de length = N.
% U es un array I \times J \times N; guarda en cada hoja un proceso 
% bidimensional definido en el rect\'angulo [a b c d],
% registrado en los tiempos de ts.

% Por defecto [a b c d]= [1 I 1 J].

% espera es el tiempo entre surf  y surf
% Por defecto espera=1/24

% Para cambiar otros par\'ametros es mejor editar y cambiarlos
% directamente al comienzo. Si se dejan como variables necesitamos
% ya una funci\'on con manual de empleo.
% 
esperadefault=1/24;
espera1=2;
toxlabel='Eje X';
toylabel='Eje Y';
totitle='Solucion en t =';
toview=[1 1 1];
toalpha=1;

if (nargin==2)
    [U,ts]=deal(varargin{:});
    espera=esperadefault;
    s=size(U);
    a=0;
    b=s(2);
    c=0;
    d=s(1);
elseif (nargin==3)
    [U,ts,espera]=deal(varargin{:});
    s=size(U);
    a=0;
    b=s(2);
    c=0;
    d=s(1);
elseif (nargin==6)
    [U,a,b,c,d,ts]=deal(varargin{:});
    espera =esperadefault;
else
    [U,a,b,c,d,ts,espera]=deal(varargin{:});
end

minU=min(U(:));
maxU=1.1*max(U(:));
if (maxU< 0)
    maxU=0;
    minU=1.1*minU;
end
    
clf
figure(1)
%set(gca,'xlim',[a b],'ylim',[c b],'zlim',[minU maxU])
pbaspect([b-a,d-c,maxU-minU])
axis([a, b, c, d, minU, maxU])
plot3([a a b b],[c c d d],[minU maxU 0 0],'.w')
grid on
hold on
surf(U(:,:,1))
shading interp
alpha(toalpha)
title([totitle,num2str(ts(1))]);
xlabel(toxlabel);
ylabel(toylabel);
view(toview)
hold off
pause
    for n=1:length(ts)
        %clf    %%MODIFICADO
        %figure(1)  %%MODIFICADO
        %set(gca,'xlim',[a b],'ylim',[c b],'zlim',[minU maxU])
        pbaspect([b-a,d-c,maxU-minU])
        axis([a, b, c, d, minU, maxU])
        plot3([a a b b],[c c d d],[minU maxU 0 0],'.w')
        caxis([0 1]);  %%MODIFICADO
        grid on
        hold on
        surf(U(:,:,n))
        shading interp
        alpha(toalpha)
        title([totitle,num2str(ts(n))]);
        xlabel(toxlabel);
        ylabel(toylabel);
        view(toview)
        if (n>1)
        pause(espera)
        else
        pause(espera1)
        end
        hold off
    end
    

end


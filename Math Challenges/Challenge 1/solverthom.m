function x = solverthom(li,ud,us,y)
% Se quiere resolver Ax=y, donde A=tridiag(a,d,b).
% Primero se ejecuta [li,ud,us]=luthomas(a,d,b), que son las variables
% que empleamos en esta funci'on, junto con el dato y.
% Vamos a sobreescribir sobre x
% Primero resolvemos Lz=y
z=y;
    for k=2:length(y)
        z(k)=y(k)-li(k-1)*z(k-1);
    end
% Ahora resolvemos Ux=z;
x=z;
x(end)=z(end)/ud(end);
    for k=length(x)-1:-1:1
        x(k)= (z(k)-us(k)*x(k+1))/ud(k);
    end
end


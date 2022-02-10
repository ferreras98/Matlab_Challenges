function [li,ud,us] = luthomas(a,d,b)
% A = tridag(a,d,b), A de tamaño n por n
% La factorizaci'on A=LU de A consta de:
%    Una matriz L con unos en la diagonal y subdiagonal descrita 
%    por un vector li de longitud n-1.
% y
%    Una matriz U con un vector ud en la diagonal, de longitud n, y de
%    un vector us en la superdiagonal, de longitud n-1.
% Todos los vectores columna y del tamaño establecido.
li=zeros(size(a));
ud=d;
li(1)=a(1)/ud(1);
ud(2)=ud(2)-li(1)*b(1);
    for k=2:length(d)-1
        li(k)=a(k)/ud(k);
        ud(k+1)=ud(k+1)-li(k)*b(k);
    end
us=b;
end
    



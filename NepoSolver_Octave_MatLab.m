%Gabriel Uri Brito Davansso
% Operacoes com numeros imaginarios
d2r = @(x) (x*pi/180);
r2d = @(x) (x*180/pi);
cis = @(x) exp(j*x);
cisd = @(x) cis(d2r(x));
angled = @(x) r2d(angle(x));

% Para declarar um numero imaginario em forma angular:
% z = r * cisd(angulo em graus)
% Para conseguir o modulo e angulo:
% [ abs(z) angled(z) ]

% ================== P1 ==================
% operador a
a = 1*cisd(120);

% matrizes A e A inversa
A = [ 1 1 1 ; 1 a**2 a; 1 a a**2; ];
Ainv = inv(A);

% ================== P2 ==================
% Algoritmo de reducao de Kron
function [Ynew] = kron_red (Y, p)
  Ynew = zeros(length(Y));
  i = 1;
  j = 1;
  for j=1:length(Y)
      for k=1:length(Y)
          if ( (j!=p) && (k!=p))
            Ynew(j,k) = Y(j,k)-(Y(j,p)*Y(p,k)/Y(p,p));
          endif
      end
  end
  Ynew([p],:)=[];
  Ynew(:,[p])=[];
endfunction

% Matriz da lista de exercicio
Y = [-5.5 1/0.4 1/0.5 0 0;
      1/0.4 -11.5 1/0.25 0 1/0.2;
      1/0.5 1/0.25 -14 1/0.125 0;
      0 0 1/0.125 -10 1/0.5;
      0 1/0.2 0 1/0.5 -7.8];
Y = i*Y
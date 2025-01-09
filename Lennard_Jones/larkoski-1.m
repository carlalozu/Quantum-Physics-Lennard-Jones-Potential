clc; clear all

% for 
n = 1
alpha = 0.5
%     if n<4
%         alpha = 0.5;
%     elseif n<7
%         alpha =0.3;
%     elseif n<8
%         alpha = 0.75;
%     else
%         alpha = 0.05;
%     end
    
mu = (12 + 24*1836)/2;
a = 4.96e7; b = 624;
l = 0; %número cuántico angular
% n = 5; %número cúantico radial

V = @(R) (l+1./2).^2./(2*mu.*R.^2)+a./(R.^12)-b./(R.^6);
R = linspace(6,10,10000);

% Derivada centrada 1
hr=R(2)-R(1);
i1=(R(2):hr:R(end));   %i+1
i3=i1-2*hr;            %i-1
dV=(V(i1)-V(i3))./(2*hr);
% Encontrar un valor de E inicial adecuado
K = [1:length(i1)-1];
E = alpha*V(i1(find(dV(K).*dV(K+1)<0)));

figure(1); plot(R,V(R),R,ones(size(R))*E,'.r')
title('Potencial Efectivo'); xlabel('Radio (bohrs)'); ylabel('Potencial')
ylim([-0.003 0.006])

% Encontar valores extramos de r tal que el radical es cero
x1 = R(find((E-V(R)).*(E-V(R+hr))<0)); %límite inferior, rango: [x1,x2]
x2 = R(find((E-V(R)).*(E-V(R+hr))<0)+1); %límite superior
xr = (x1+x2)/2;

% Regla del Simpson 1/3 (múltiples segmentos, números impares)
EE = linspace(E*2,0,1000);
k = linspace(xr(1),xr(2),1000);
%k2 = linspace(xr(1),xr(2),20000);
xi = k(2:2:(end-1));
xj = k(3:2:(end-2));

for t = 1:length(EE)
    % Integral
    %It(t) = (k(2)-k(1)).*(sqrt(f(xr(1),EE(t)))+4*sum(sqrt(f(xi,EE(t))))+2*sum(sqrt(f(xj,EE(t))))+sqrt(f(xr(2),EE(t))))./3;
    It(t) = (k(2)-k(1)).*(sqrt(EE(t)-V(xr(1)))+4*sum(sqrt(EE(t)-V(xi)))+2*sum(sqrt(EE(t)-V(xj)))+sqrt(EE(t)-V(xr(2))))./3;
    
    % Trapz
    %It(t) = trapz(k,sqrt(EE(t)-V(k)));
    %Is2(t) = trapz(k2,sqrt(f(k2,EE(t))));
    % Romberg 
    %It(t) = 4/3*Is2(t)-1/3*Is(t);
    % Expresión completa
    p(t) = sqrt(2*mu)./pi.*It(t)-(n+1/2);
end
%
figure(2); plot(EE,real(p))
title('f(E)'); xlabel('Energía'); ylabel('Valor de f(E)')
k = 1:length(p)-1;
Ef = EE(find(real(p(k)).*real(p(k+1))<0));
% end

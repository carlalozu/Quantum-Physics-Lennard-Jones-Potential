clear all

mu = (12 + 24*1836)/2;
a = 4.96e7; b = 624;
l = 1; %número cuántico angular

% Potencial
V = @(R) (l.^2+l+1./4)./(2*mu.*R.^2)+a./(R.^12)-b./(R.^6);

% Valor de E inicial
R = linspace(6,20,100000);
hr=R(2)-R(1);
i1=(R(2):hr:R(end));   %i+1
dV=(V(i1)-V(i1-2*hr))./(2*hr);
K = 1:length(i1)-1;
E = 0.5*V(i1(dV(K).*dV(K+1)<0));

figure(2); plot(R,V(R),R,ones(size(R))*E,'.r')
title('Potencial Efectivo'); xlabel('Radio (bohrs)'); ylabel('Potencial')
ylim([-0.0021 0])

% Integración trapezoidal
EE = linspace(2*E,-0.00009,1000);
It = zeros(1,length(EE));
f = It;
Ef = zeros(1,11);

Ef = [-0.001795582992563
  -0.001488869064956
  -0.001217319791470
  -0.000978981580209
  -0.000773854431173
  -0.000598031160570
  -0.000451511768402
  -0.000330389070876
  -0.000230755884201
  -0.000154565800273
  -0.000097911635301].';

%%

for n = 0:10 %número cúantico radial
    for t = 1:length(EE)
        % Encontar valores extramos de r tal que el radical es cero
        xr = R((EE(t)-V(R)).*(EE(t)-V(R+hr))<0);

        k=linspace(xr(1),xr(2),1000);
        % Integral y función completa
        It(t) = trapz(k,sqrt(EE(t)-V(k)));
        f(t) = sqrt(2*mu)/pi*It(t)-(n+1/2);
    end

    k = 1:length(f)-1;
    Ef(n+1) = EE(real(f(k)).*real(f(k+1))<0);
end
%%

figure(1); plot(R,V(R),'k')
hold on

for n=0:10
NivelE = NaN(1,size(R,2));
NivelE(V(R) <= Ef(n+1)) =  Ef(n+1);
p(n+1) = plot(R,NivelE);
end

ylim([-0.0021 0]); 
grid on
title('Niveles de energía para la molécula de Mg2 para el potencial 6-12'); xlabel('Radio (Bohrs)'); ylabel('Potencial (Hartree)');
legend([p(1) p(2) p(3) p(4) p(5) p(6) p(7) p(8) p(9) p(10) p(11)],{'n=0','n=1','n=2','n=3','n=4','n=5','n=6','n=7','n=8','n=9','n=10'},'Location','east');


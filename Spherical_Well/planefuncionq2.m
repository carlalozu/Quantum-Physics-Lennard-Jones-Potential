clc; clear all; close all;
% condiciones inciales
a = 1; h = 1; M = 1/2; 

% mandamos a llamar los coeficientes
S = load('coeff2');
coeff2 = cell2mat(struct2cell(S));

% generación del espacio
N = 2^6;
rho = linspace(0,a,N);
az = linspace(0,2*pi,N);
theta = pi/2;
[phi,r] = meshgrid(az,rho);

% descomposición en eigenestados
L = 5; nmax = 5; qt=0;

k= 50; seg = 20;
saveq = zeros(k,N,N);
p = 1;
for t = 0:seg/k:seg %tiempo
for n = 1:nmax
for l=0:L
    
    i=1;
    % bessel esférica
    J = sqrt(pi./(2.*betaln(l,n).*r)).*besselj(l+1/2,betaln(l,n).*r);
    % normalización
    Aln = sqrt(2)/(sqrt(pi./(2.*betaln(l,n))).*besselj(l+3/2,betaln(l,n)));
    % energía
    En = betaln(l,n).^2;
    
    for m = -l:l
    % eigenstados
    Psi = Aln.*J.*ylm2d(l,m,phi,theta);
    qt = qt + coeff2(n,i,l+1).*Psi.*exp(-1i*En*t);
    
    i=i+1;
    end
end
end
saveq(p,:,:) = qt;
disp(p)
p = p + 1;
end

%% Gráfica

X = r.*cos(phi);
Y = r.*sin(phi);

for p = 1:k+1
figure(1);
kj = reshape(saveq(p,:,:),64,64);
surf(X,Y,abs(kj).^2);
view(2)
set(gca,'FontSize',12);
set(0,'defaultTextInterpreter','latex');
title(['$|Q(t=$' num2str(seg/k*(p-1)) '$)|^2$'])
axis equal
xlabel('X'); ylabel('Y');
colormap parula
shading interp
pause(0.1)
end

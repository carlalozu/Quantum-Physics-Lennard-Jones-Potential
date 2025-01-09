clc; clear all; close all;
% condiciones inciales
a = 1; h = 1; M = 1/2; 

% mandamos a llamar los coeficientes del problema pasado
S = load('coeff1');
coeff1 = cell2mat(struct2cell(S));

% generación del espacio
N=2^12;
r = linspace(0.0000001,a,N);
B = 0.32792284570150765;
Q = B.*(1-r.^2);

% decomposición in eigenestados
L=15; nmax = 10;
coeff2= zeros(nmax,size(-L:L,2),L+1); qY = 0;

for n = 1:nmax
    for l=0:L
    i=1;
    for m = -l:l
    % bessel esférica
    J = sqrt(pi./(2.*betaln(l,n).*r./a)).*besselj(l+1/2,betaln(l,n).*r./a);
    % normalización
    Aln = sqrt(2/a^3)/(sqrt(pi./(2.*betaln(l,n))).*besselj(l+1+1/2,betaln(l,n)));
    % eigenstates
    Psi = Aln.*J;

    coeff2(n,i,l+1) = trapz(conj(Psi).*Q.*r.^2.*coeff1(i,l+1)).*(a./N);
    qY = qY + coeff2(n,i,l+1).*Psi;
    disp(coeff2(n,i,l+1))
    i=i+1;
    end
    end
end

% reshape(coeff2(n,:,:),size(-L:L,2),L+1)   para visualizar los coeffs

% el coeficiente máximo corresponde a 
MaxC = max(max(max(coeff2)));
% que es n=1,l=0,m=0

%%
Et = 0;
for n = 1:nmax
    for l=0:L
        i=1;
    for m = -l:l
        %para el valor esperado de la energía
        En = h^2/(2*M*a)*betaln(l,n)^2.*conj(coeff2(n,i,l+1)).*coeff2(n,i,l+1);
        Et = Et + En;

        %para el valor esperado del momentum angular L^2
        Ln = l*(l+1).*conj(coeff2(n,i,l+1)).*coeff2(n,i,l+1);
        Lt = Lt + Ln;

        i=i+1;
        end
    end
end
function Ylm = ylm(l,M,phi,theta)
m = abs(M);
L = legendre(l,cos(theta));
if l == 0 
    Lm = L;
else
    Lm = reshape(L(m+1,:,:),size(phi));
end
norm = sqrt((2*l+1).*factorial(l-m)/(4*pi*factorial(l+m)));
phase = (-1)^m.*exp(1i.*m.*phi);

Ylm = norm.*Lm.*phase;
if M<0
    Ylm = (-1)^m.*conj(Ylm);
end
end

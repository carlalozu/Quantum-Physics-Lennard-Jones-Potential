function betal = betaln(l,n)
N=2^17; up = 100; h = up/N;
r = linspace(0,up,N);
j = @(l,x) sqrt(pi./(2.*x)).*besselj(l+1/2,x);
betal = r(j(l,r).*j(l,r+h)<0) + h/2;
betal = betal(n);
end

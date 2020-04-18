%Cette fonction et presque la meme que lnf1 la seule difference ceest que 
%j ai ajoute un autre variable, x=xa afin de tracer la courbe x=f(DAB)    
function [dt]= lnf2 ( x,s1,s2)
d21=2.10*10^-5;
d12=2.67*10^-5;
x(2,1)=x;
x(1,1)=1-x(2,1);
q(1,1)=1.4;
q(2,1)=1.432;
a(2,1)=s1;
a(3,1)=s2;
a(1,1)=0;
a(4,1)=0;
T=313.13;
r(1,1)=1.4311;
r(2,1)=0.92;
l=0;
m=0;
for z=1:1:2
         l(z,1)=l+(x(z,1).*q(z,1));
         y(z,1)=(r(z,1)).^(1/3);
         m(z,1)=m+(x(z,1).*y(z,1));
         for v=1:1:2
      p(v,1)=(x(v,1).*q(v,1))./l(z,1);
      f(v,1)=(x(v,1).*r(v,1).^(1/3))/m(z,1);
         end
end
for g=1:1:4
    t(g,1)=exp(-a(g,1)./T);
end

 h(1,1)=p(1,1)*t(2,1)+ p(2,1)*t(1,1);
 h(2,1)=p(2,1)*t(4,1)+ p(1,1)*t(3,1);
 
 c(1,1)=p(1,1)*t(1,1)/h(1,1);
 c(2,1)=p(2,1)*t(2,1)/h(1,1);
 c(3,1)=p(2,1)*t(4,1)/h(2,1);
 c(4,1)=p(1,1)*t(3,1)/h(2,1);
 
lnDT=x(1,1)*log(d21)+x(2,1)*log(d12)+2*sum((x.*log(x./f)))+2*prod(x)*sum((f./x).*(1-y./flip(y)))+x(1,1)*q(2,1)*((1-c(4,1)^2)*log(t(3,1))+(1-c(1,1)^2)*t(2,1)*log(t(2,1)))+x(2,1)*q(1,1)*((1-c(2,1)^2)*log(t(2,1))+(1-c(3,1)^2)*t(3,1)*log(t(3,1)));
dt=exp(lnDT);
end
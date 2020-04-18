%C est la fonction qui va me permettre de déterminer a12 et a21
%C est une fonction a deux variables  %a12=s1%, %a21=s2%
function [lnFC]= lnf1 ( s1 , s2 )
% 1/Cette premiere partie concerne l initialisation des donnes
%Puisque j ai  travaille avec les vecteurs vous trouviez que :
%x1=xb=x(1,2)
%x2=xa=x(2,1)
%qb=q1=q(1,1)
%qa=q2=q(2,1)
%a(2,1)=aAB (a21)
%a(3,1)=aBA (a12)
%a(1,1)=a11 a(BB)
%a(4,1)=a22 a(AA)
%r2=ra=r(1,1)
%r1=rb=r(2,1)
d21=2.10*10^-5;
d12=2.67*10^-5;
x(1,1)=0.65;
x(2,1)=0.35;
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
%Afin de calcule i et  j ai utilise  les boucles 
% avec p(v,1)= phi(i)
%f(v,1)= teta(j)
%y(z,1)= londa(i)                          
for z=1:1:2
         l(z,1)=l+(x(z,1).*q(z,1));
         y(z,1)=(r(z,1)).^(1/3);
         m(z,1)=m+(x(z,1).*y(z,1));
         for v=1:1:2
      p(v,1)=(x(v,1).*q(v,1))./l(z,1);
      f(v,1)=(x(v,1).*r(v,1).^(1/3))/m(z,1);
         end
end
% tau(ji)=t(g,1)
for g=1:1:4
    t(g,1)=exp(-a(g,1)./T);
end
% teta(ji)
 h(1,1)=p(1,1)*t(2,1)+ p(2,1)*t(1,1);
 h(2,1)=p(2,1)*t(4,1)+ p(1,1)*t(3,1);
 
 c(1,1)=p(1,1)*t(1,1)/h(1,1);
 c(2,1)=p(2,1)*t(2,1)/h(1,1);
 c(3,1)=p(2,1)*t(4,1)/h(2,1);
 c(4,1)=p(1,1)*t(3,1)/h(2,1);
 
lnFC=x(1,1)*log(d21)+x(2,1)*log(d12)+2*sum((x.*log(x./f)))+2*prod(x)*sum((f./x).*(1-y./flip(y)))+x(1,1)*q(2,1)*((1-c(4,1)^2)*log(t(3,1))+(1-c(1,1)^2)*t(2,1)*log(t(2,1)))+x(2,1)*q(1,1)*((1-c(2,1)^2)*log(t(2,1))+(1-c(3,1)^2)*t(3,1)*log(t(3,1)));
end
clc
clear
D21=1.33*10^-5;
 %Puis s1 et s2 est une chaîne de caracteres, on utilisant  syms x et y vont
%devenir  des variables symboliques   
syms s1 s2 
var=[s1,s2];
%la fonction Matlab utilisee pour resoudre une equation transcendante
%contient une ou plusieurs fonctions  transcendante comme
%(sin(x), log(x) ou ex , c’est la fonction solve.
sol=vpasolve((abs(lnf1(s1,s2)-log(D21)))==0,var);
%Puisque tout le calcule  est symbolique pour avoir les valeurs numerique 
%de a12 et a21 j ai utilisé les fonctions double ainsi num2str(Numero to string ) 
disp(['s1= aAB (a21) = ',num2str(double(sol.s1))])
disp(['s2 =aBA (a12)= ',num2str(double(sol.s2))])
disp(['err = ',num2str(abs(lnf1(double(sol.s1),double(sol.s2))-log(D21))/log(D21))])
disp(['DAB= ',num2str(exp(lnf1(double(sol.s1),double(sol.s2))))])
%C est la boucle qui va me permettre de tracer la courbe  DAB=f(Xa)
x=0:0.01:0.7;
n=0.7/0.01;
for i=1:n+1
    y(i)=lnf2(x(i),double(sol.s1),double(sol.s2));
end
plot(x,y,'-o','MarkerIndices',1:2:length(y))
ylabel('DAB')
xlabel('xa')
%Puisque  la valeur de l erreur égale à 0  donc les valeurs des a12 a21 son
%juste.
%l ecart entre la valeur théorique et experimentatl de mon programme egale
%a 0.

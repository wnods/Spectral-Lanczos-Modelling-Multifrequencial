%======================================================================
%======================================================================
%Plotagem dos campos estimados para uma unica frequencias base,
%======================================================================
%======================================================================
close all;
clear all;
clc;

%Entrada de dados

ind   = load('indices_perfil.dat'); 
coord = load('coordenadas.dat');

x = coord(:,1);
nx=length(ind);
%=======================================
%Dados calculados freq 100
%=======================================


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%50hz%%%%%%%%%%%%%%%
IMAG100 = load('IMAG100.dat');
REAL100 = load('REAL100.dat');

campo_Cal100 = load('Campo_100.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);



vec1_100 =  REAL100+1i*IMAG100;
vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(1) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
plot(x,REAL100(ind),'b',x,campo_Cal100(ind,1),'-ob',...
     x,IMAG100(ind),'r',x,campo_Cal100(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
campo_Cal100 = load('Campo_100_151.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);


vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(2) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
campo_Cal100 = load('Campo_100_252.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);


vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(3) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campo_Cal100 = load('Campo_100_302.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);




vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(4) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
campo_Cal100 = load('Campo_100_504.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);


vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(5) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
campo_Cal100 = load('Campo_100_655.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);


vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(6) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
campo_Cal100 = load('Campo_100_756.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);


vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(7) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
campo_Cal100 = load('Campo_100_1260.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);


vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(8) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

campo_Cal100 = load('Campo_100_2520.dat');

field2 = campo_Cal100(ind,1)+sqrt(-1)*campo_Cal100(ind,2);


vec2_100 = field2;


%Calculo do erro

k1=1;
vec = 0;
for j2=1:nx

n1 = sqrt(real(vec1_100(ind(j2)))^2 + imag(vec1_100(ind(j2)))^2) ;
n2 = sqrt(real(vec2_100((j2)))^2 + imag(vec2_100((j2)))^2) ;
 
vec(k1) =  (n1-n2);

k1=k1+1;
end 

err(9) = sqrt(dot(vec,vec)/(k1-1)); %Erro Quadratico




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x_1=[1:9];

err1(1) = err(1);
err1(2) = err(2);
err1(3) = err(4);
err1(4) = err(6);
err1(5) = err(8);


x_2=[1:5]

figure(2)
plot(x_1,err,'*-b')
legend('erro')
set(gca,'fontsize',16)
title('Est com pod1');
xlabel('x (m)')
grid;




figure(21)
plot(x_2,err1,'*-b')
legend('erro')
set(gca,'fontsize',16)
title('Est com pod1');
xlabel('x (m)')
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(x,REAL100(ind),'b',x,campo_Cal100(ind,1),'-ob',...
     x,IMAG100(ind),'r',x,campo_Cal100(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;

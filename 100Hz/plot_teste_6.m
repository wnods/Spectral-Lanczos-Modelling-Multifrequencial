%======================================================================
%======================================================================
%Plotagem dos campos estimados para uma unica frequencias base,
%mas usando ponderador
%Figuras salvas no formato em eps e jpg.
%Modelo de dimensões 5041x5041
%Frequencias variando de 50 a 150Hz
%Erro calculado e cada campo salvo por vetores diferentes
%======================================================================
%======================================================================
close all;
clear all;
clc;

%Entrada de dados


ind   = load('indices_perfil.dat'); 
coord = load('coordenadas.dat');



load imagB100.dat;  %vetor fonte: part imag da freq base
load realB100.dat; %vetor fonte: part real da freq base



x = coord(:,1);
nx=length(ind);

realB_fb =realB100 ;
imagB_fb =imagB100;

B_fb = realB_fb +sqrt(-1)*imagB_fb;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%50hz%%%%%%%%%%%%%%%
IMAG50 = load('IMAG50.dat');
REAL50 = load('REAL50.dat');

campo_Cal50 = load('Campo_50_655.dat');

load imagB50.dat; %vetor fonte: part imag da outra freq
load realB50.dat;%vetor fonte: part real da outra freq

realB_f=realB50;
imagB_f=imagB50;

B_f = realB_f +sqrt(-1)*imagB_f;


mat1 = B_f*transp(B_f);

mat2 = pinv(mat1);

pond1 = transp(B_f)*mat2*B_fb; 

field2 = campo_Cal50(ind,1)+sqrt(-1)*campo_Cal50(ind,2);

field3 = (1/(pond1))*field2;

vec1 =  REAL50+1i*IMAG50;
vec2 = field3;
vec3 = field2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(50)
plot(x,REAL50(ind),'b',x,campo_Cal50(ind,1),'-ob',...
     x,IMAG50(ind),'-.r',x,campo_Cal50(ind,2),'-*r')
set(gca,'fontsize',16)
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('X [m]')
ylabel('Amplitude [a. u.]')
grid;

figure(501)
plot(x,REAL50(ind),'b',x,real(field3),'-ob',...
     x,IMAG50(ind),'-.r',x,imag(field3),'-*r')
set(gca,'fontsize',16)
legend('real Exact','real Cal Pond','imag Exact','imag Cal Pond')
xlabel('X [m]')
ylabel('Amplitude [a. u.]')
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%60hz%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%50hz%%%%%%%%%%%%%%%
IMAG150 = load('IMAG150.dat');
REAL150 = load('REAL150.dat');

campo_Cal150 = load('Campo_150_655.dat');

load imagB150.dat; %vetor fonte: part imag da outra freq
load realB150.dat;%vetor fonte: part real da outra freq

realB_f=realB150;
imagB_f=imagB150;

B_f = realB_f +sqrt(-1)*imagB_f;


mat1 = B_f*transp(B_f);

mat2 = pinv(mat1);

pond1 = transp(B_f)*mat2*B_fb; 

field2 = campo_Cal150(ind,1)+sqrt(-1)*campo_Cal150(ind,2);

field3 = (1/(pond1))*field2;

vec1 =  REAL150+1i*IMAG150;
vec2 = field3;
vec3 = field2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(150)
plot(x,REAL150(ind),'b',x,campo_Cal150(ind,1),'-ob',...
     x,IMAG150(ind),'-.r',x,campo_Cal150(ind,2),'-*r')
set(gca,'fontsize',16)
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('X [m]')
ylabel('Amplitude [a. u.]')
grid;

figure(1501)
plot(x,REAL150(ind),'b',x,real(field3),'-ob',...
     x,IMAG150(ind),'-.r',x,imag(field3),'-*r')
set(gca,'fontsize',16)
legend('real Exact','real Cal Pond','imag Exact','imag Cal Pond')
xlabel('X [m]')
ylabel('Amplitude [a. u.]')
grid;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%60hz%%%%%%%%%%%%%%%



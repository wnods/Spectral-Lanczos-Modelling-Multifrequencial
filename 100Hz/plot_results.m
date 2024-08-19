close all;
clear all;
ind   = load('indices_perfil.dat'); 
coord = load('coordenadas.dat');
%=======================================
%Dados para frequencia base de 80hz
%=======================================
for i=1:11
freq(i) = 50 +(i-1)*10;
end

x = coord(:,1);
nx=length(ind);


load imagB100.dat;  %vetor fonte: part imag da freq base
load realB100.dat; %vetor fonte: part real da freq base

realB_fb =realB100;
imagB_fb =imagB100;

B_fb = realB_fb +sqrt(-1)*imagB_fb;


IMAG50 = load('IMAG50.dat');
REAL50 = load('REAL50.dat');
campo_Cal50 = load('Camp_est_1.dat');
load imagB50.dat; 
load realB50.dat;



IMAG60 = load('IMAG60.dat');
REAL60 = load('REAL60.dat');
campo_Cal60 = load('Camp_est_2.dat');
load imagB60.dat; 
load realB60.dat;

IMAG70 = load('IMAG70.dat');
REAL70 = load('REAL70.dat');
campo_Cal70 = load('Camp_est_3.dat');
load imagB70.dat; 
load realB70.dat;

IMAG80 = load('IMAG80.dat');
REAL80 = load('REAL80.dat');
campo_Cal80 = load('Camp_est_4.dat');
load imagB80.dat; 
load realB80.dat;


IMAG90 = load('IMAG90.dat');
REAL90 = load('REAL90.dat');
campo_Cal90 = load('Camp_est_5.dat');
load imagB90.dat; 
load realB90.dat;


IMAG100 = load('IMAG100.dat');
REAL100 = load('REAL100.dat');
campo_Cal100 = load('Camp_est_6.dat');


IMAG110 = load('IMAG110.dat');
REAL110 = load('REAL110.dat');
campo_Cal110 = load('Camp_est_7.dat');
load imagB110.dat; 
load realB110.dat;

IMAG120 = load('IMAG120.dat');
REAL120 = load('REAL120.dat');
campo_Cal120 = load('Camp_est_8.dat');
load imagB120.dat; 
load realB120.dat;

IMAG130 = load('IMAG130.dat');
REAL130 = load('REAL130.dat');
campo_Cal130 = load('Camp_est_9.dat');
load imagB130.dat; 
load realB130.dat;

IMAG140 = load('IMAG140.dat');
REAL140 = load('REAL140.dat');
campo_Cal140 = load('Camp_est_10.dat');
load imagB140.dat; 
load realB140.dat;


IMAG150 = load('IMAG150.dat');
REAL150 = load('REAL150.dat');
campo_Cal150 = load('Camp_est_11.dat');
load imagB150.dat; 
load realB150.dat;


x = coord(:,1);

nx=length(ind);


realB_f=realB50;
imagB_f=imagB50;

B_f = realB_f +sqrt(-1)*imagB_f;


mat1 = B_f*transp(B_f);
mat2 = pinv(mat1);
pond1 = transp(B_f)*mat2*B_fb; 

field1 = campo_Cal50(ind,1)+sqrt(-1)*campo_Cal50(ind,2);
field2 = (1/(pond1))*field1;






vec1 =  REAL60(ind)+1i*IMAG60(ind);
vec2 = campo_Cal60(ind,1)+1i*campo_Cal60(ind,2);

Abs_ext(2) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(2) = sqrt((dot(vec2,vec2))/(nx-1));


vec1 =  REAL70(ind)+1i*IMAG70(ind);
vec2 = campo_Cal70(ind,1)+1i*campo_Cal70(ind,2);

Abs_ext(3) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(3) = sqrt((dot(vec2,vec2))/(nx-1));

vec1 =  REAL80(ind)+1i*IMAG80(ind);
vec2 = campo_Cal80(ind,1)+1i*campo_Cal80(ind,2);

Abs_ext(4) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(4) = sqrt((dot(vec2,vec2))/(nx-1));

vec1 =  REAL90(ind)+1i*IMAG90(ind);
vec2 = campo_Cal90(ind,1)+1i*campo_Cal90(ind,2);

Abs_ext(5) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(5) = sqrt((dot(vec2,vec2))/(nx-1));

vec1 =  REAL100(ind)+1i*IMAG100(ind);
vec2 = campo_Cal100(ind,1)+1i*campo_Cal100(ind,2);

Abs_ext(6) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(6) = sqrt((dot(vec2,vec2))/(nx-1));


vec1 =  REAL110(ind)+1i*IMAG110(ind);
vec2 = campo_Cal110(ind,1)+1i*campo_Cal110(ind,2);

Abs_ext(7) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(7) = sqrt((dot(vec2,vec2))/(nx-1));

vec1 =  REAL120(ind)+1i*IMAG120(ind);
vec2 = campo_Cal120(ind,1)+1i*campo_Cal120(ind,2);

Abs_ext(8) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(8) = sqrt((dot(vec2,vec2))/(nx-1));


vec1 =  REAL130(ind)+1i*IMAG130(ind);
vec2 = campo_Cal130(ind,1)+1i*campo_Cal130(ind,2);

Abs_ext(9) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(9) = sqrt((dot(vec2,vec2))/(nx-1));



vec1 =  REAL140(ind)+1i*IMAG140(ind);
vec2 = campo_Cal140(ind,1)+1i*campo_Cal140(ind,2);

Abs_ext(10) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(10) = sqrt((dot(vec2,vec2))/(nx-1));



vec1 =  REAL150(ind)+1i*IMAG150(ind);
vec2 = campo_Cal150(ind,1)+1i*campo_Cal150(ind,2);

Abs_ext(11) = sqrt((dot(vec1,vec1))/(nx-1));
Abs_cal(11) = sqrt((dot(vec2,vec2))/(nx-1));



x1=0;
for i=1:11
freq(i) = 50 + (i-1)*x1;
x1=10;
end


err_freq = 100*abs(abs(Abs_ext)-abs(Abs_cal))./abs(Abs_ext);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(50)
plot(x,REAL50(ind),'b',x,campo_Cal50(ind,1),'-ob',...
     x,IMAG50(ind),'r',x,campo_Cal50(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;

figure(501)
plot(x,REAL50(ind),'b',x,real(field2),'-ob',...
     x,IMAG50(ind),'r',x,imag(field2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;



figure(60)
plot(x,REAL60(ind),'b',x,campo_Cal60(ind,1),'-ob',...
     x,IMAG60(ind),'r',x,campo_Cal60(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;

figure(70)
plot(x,REAL70(ind),'b',x,campo_Cal70(ind,1),'-ob',...
     x,IMAG70(ind),'r',x,campo_Cal70(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;

figure(80)
plot(x,REAL80(ind),'b',x,campo_Cal80(ind,1),'-ob',...
     x,IMAG80(ind),'r',x,campo_Cal80(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;

figure(90)
plot(x,REAL90(ind),'b',x,campo_Cal90(ind,1),'-ob',...
     x,IMAG90(ind),'r',x,campo_Cal90(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;

figure(100)
plot(x,REAL100(ind),'b',x,campo_Cal100(ind,1),'-ob',...
     x,IMAG100(ind),'r',x,campo_Cal100(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;
hold on

figure(110)
plot(x,REAL110(ind),'b',x,campo_Cal110(ind,1),'-ob',...
     x,IMAG110(ind),'r',x,campo_Cal110(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;
hold on


figure(120)
plot(x,REAL120(ind),'b',x,campo_Cal120(ind,1),'-ob',...
     x,IMAG120(ind),'r',x,campo_Cal120(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;
hold on

figure(130)
plot(x,REAL130(ind),'b',x,campo_Cal130(ind,1),'-ob',...
     x,IMAG130(ind),'r',x,campo_Cal130(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;
hold on

figure(140)
plot(x,REAL140(ind),'b',x,campo_Cal140(ind,1),'-ob',...
     x,IMAG140(ind),'r',x,campo_Cal140(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;
hold on

figure(150)
plot(x,REAL150(ind),'b',x,campo_Cal150(ind,1),'-ob',...
     x,IMAG150(ind),'r',x,campo_Cal150(ind,2),'-or')
legend('real Exact','real Cal','imag Exact','imag Cal')
xlabel('x (m)')
grid;
hold on



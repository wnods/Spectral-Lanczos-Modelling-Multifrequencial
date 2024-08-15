clear all

Hy  = load('Hyfield.dat');
Exz = load('ExEzfield.dat');
res = load('Res_apar.dat');

N = length(Hy(:,1));

x = zeros(N,1);

xi = -8000;
xf =  8000;
dx = (xf - xi)/(N-1);

for i = 1:N
   x(i) = xi + (i-1)*dx; 
end

lw = 1;
fs = 13;

figure(1)
subplot(2,1,1)
plot(x,Hy(:,1),'-o',x,Hy(:,2),'-o','Linewidth',lw)
set(gca,'Fontsize',fs)
xlabel('x (m)')
ylabel('H_y^s (A/m)')
xlim([-3000 3000])
legend('Re','Im')
hold on
grid on
%grid

%figure
subplot(2,1,2)
plot(x,Exz(:,1),'-o',x,Exz(:,2),'-o','Linewidth',lw)
set(gca,'Fontsize',fs)
xlabel('x (m)')
ylabel('E_x^s (V/m)')
xlim([-3000 3000])
legend('Re','Im')
hold on
grid on
%grid

% figure
% subplot(3,1,3)
% plot(x,Exz(:,3),'-or',x,Exz(:,4),'-ob','Linewidth',lw)
% set(gca,'Fontsize',fs)
% xlabel('x (m)')
% ylabel('E_z^s (V/m)')
% legend('Re','Im')
% grid

% figure
% subplot(2,1,1)
% plot(x,res(:,1),'-or','Linewidth',lw)
% set(gca,'Fontsize',fs)
% xlabel('x (m)')
% ylabel('\rho_a (\Omega m)')
% legend('Resistividade')
% grid
% subplot(2,1,2)
% plot(x,res(:,2),'-or','Linewidth',lw)
% set(gca,'Fontsize',fs)
% xlabel('x (m)')
% ylabel('\phi (Graus)')
% legend('Fase')
% grid


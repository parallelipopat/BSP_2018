T = 20;
N = 1000;
h = T/N;

Z_exp = zeros(2,N); Z_imp = zeros(2,N);Z_symp = zeros(2,N);
f = @(x) [x(1)-x(1)*x(2); x(1)*x(2)-2*x(2)];
df = @(x, y) [1-x(2); x(1) - 2];

F_IE = @(x, y0, h) x - y0 - h*f(x);
DF_IE = @(x, h) eye(2) - h*df(x);
N_IE = @(y0, h) y0+h - (DF_IE(y0+h, h))^-1*F_IE(y0+h,y0,h);

ysolv_EE = [2;4]; ysolv_IE = [2;4]; ysolv_symp = [2;4];
Z_exp(:,1) = ysolv_EE; Z_imp(:,1) = ysolv_IE; Z_symp(:,1) = ysolv_symp;

for ii = 2: N
   ysolv_EE = ysolv_EE + h*f(ysolv_EE);
   ysolv_IE = N_IE(ysolv_IE, h);
   
   Z_symp(1,ii) = Z_symp(1,ii-1)/(1-h*(1-Z_symp(2,ii-1)));
   Z_symp(2,ii) = Z_symp(2,ii-1)*(1+h*(Z_symp(1,ii)-2));

   Z_exp(:,ii) = ysolv_EE;
   Z_imp(:,ii) = ysolv_IE;
   
end

x = 0.2:0.2:max(Z_exp(1,:));
y = 0.2:0.2:max(Z_exp(2,:));
[X,Y] = meshgrid(x,y);
Z = X - 2*real(log(X)) + Y - real(log(Y));

ax1 = subplot(1,3,1);
plot(ax1, Z_exp(1,:), Z_exp(2,:),'b')
hold on;
plot(ax1, [2], [4], 'o');
hold on;
contour(ax1, X,Y,Z,'k:');
title(ax1,'Explicit Euler');
xlabel(ax1, 'u');
ylabel(ax1, 'v');

ax2 = subplot(1,3,2);
plot(ax2, Z_imp(1,:), Z_imp(2,:),'b')
hold on;
plot(ax2, [2], [4], 'o');
hold on;
contour(ax2, X,Y,Z,'k:');
title(ax2,'Implicit Euler');
xlabel(ax2, 'u');
ylabel(ax2, 'v');

ax3 = subplot(1,3,3);
plot(ax3, Z_symp(1,:), Z_symp(2,:),'b')
hold on;
plot(ax3, [2], [4], 'o');
hold on;
contour(ax3, X,Y,Z,'k:');
title(ax3,'Symplectic Euler');
xlabel(ax3, 'u');
ylabel(ax3, 'v');

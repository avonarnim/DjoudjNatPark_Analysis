r1 = 1.4;
k1 = 50;
alpha = 0.16;
beta = 0.015;
D = 30;
h1 = @(x) 0.3+0.1*cos(2*pi*x);
h2 = @(x) 0.3+0.07*sin(2*pi*x);
h2 = @(x) 0.3+0.07*sin(2*pi*x);
Lambda = @(x) 1.2+sin(x)*exp(sin(x+pi/4))/(exp(cos(x+pi/4)));
lambda = @(x) 0.05+0.012*cos(2*pi*x);
L = @(x) 0.6+cos(x)*exp(cos(x+7*pi/4))/(exp(sin(x+7*pi/4)));
d1 = 1.5;
Cb = 30;
gamma = 0.1;

timesteps = 1000;

trackX = zeros(1,timesteps);
trackY = zeros(1,timesteps);
trackZ = zeros(1,timesteps);


Xn = 30;
Yn = 0.5;
Zn = 15;

Xplus = Xn;
Yplus = Yn;
Zplus = Zn;

for index=1:timesteps
    t = (index-1)*20/timesteps;

    trackX(index) = Xn;
    trackY(index) = Yn;
    trackZ(index) = Zn;

    XDeriv = r1*Xn*(1-Xn/k1) - alpha*Xn*Yn/(Xn+D) - h1(t)*Xn
    YDeriv = -d1*Yn + beta*Xn*Yn/(Xn+D) - h2(t)*Yn - lambda(t)*Yn*Zn + Lambda(t) - L(t)*Yn
    ZDeriv = Cb*Yn - gamma*Zn^2

    Xplus = Xn + XDeriv * 20/timesteps;
    Yplus = Yn + YDeriv * 20/timesteps;
    Zplus = Zn + ZDeriv * 20/timesteps;

    Xn = Xplus;
    Yn = Yplus;
    Zn = Zplus;
end


t = tiledlayout(3,1);

nexttile

plot(trackX, '-');
title('Prey Population')
xlabel('Timestep');
ylabel('Population');

nexttile

plot(trackY, '-');
title('Predator Population');
xlabel('Timestep');
ylabel('Population');

nexttile

plot(trackZ, '- ');
title('Tourist Population');
xlabel('Timestep');
ylabel('Population');


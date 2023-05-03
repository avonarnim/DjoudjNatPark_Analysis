r1 = 1.4;
k1 = 50;
alpha = 0.16;
beta = 0.015;
D = 30;
h1 = @(x) 0.1+0.05*cos(2*pi*x);
h2 = @(x) 0.3+0.07*sin(2*pi*x);
Lambda = @(x) 4.8+2*sin(2*pi*x);
lambda = @(x) 0.05+0.012*cos(2*pi*x);
L = @(x) 1.1+0.7*cos(2*pi*x);
d1 = 1.5;
Cb = 30;
Cf = 0.5; % assume that fish population affects tourism too. Note that when Cf=0, the tourist population is not affected by fish
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

    % Zn factor of h1(t)XnZn term imparts tourism impact on fishing
    % frequency
    XDeriv = r1*Xn*(1-Xn/k1) - alpha*Xn*Yn/(Xn+D) - h1(t)*Xn*Zn;
    YDeriv = -d1*Yn + beta*Xn*Yn/(Xn+D) - (h2(t)+lambda(t))*Yn*Zn + Lambda(t) - L(t)*Yn;
    ZDeriv = Cb*Yn + Cf*Xn - gamma*Zn^2;

    Xplus = Xn + XDeriv * 20/timesteps;
    Yplus = Yn + YDeriv * 20/timesteps;
    Zplus = Zn + ZDeriv * 20/timesteps;

    Xn = Xplus;
    Yn = Yplus;
    Zn = Zplus;
end


tiles = tiledlayout(3,1);
t = linspace(0,20,timesteps); 

nexttile

plot(t, trackX, '-');
title('Prey Population')
xlabel('Year');
ylabel('Population');

nexttile

plot(t, trackY, '-');
title('Predator Population');
xlabel('Year');
ylabel('Population');
ylim([0,1]);

nexttile

plot(t, trackZ, '- ');
title('Tourist Population');
xlabel('Year');
ylabel('Population');


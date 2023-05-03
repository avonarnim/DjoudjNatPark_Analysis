r1 = 3;
r2 = 2;
k1 = 1;
k2 = 1;
alpha = 0.16;
beta = 0.15;
D = 0.5;
h1 = @(x) 0.1+0.1*cos(2*pi*x);
h2 = @(x) 0.07+0.07*sin(2*pi*x);
lambda = @(x) 0.012+0.012*cos(2*pi*x);
Lambda = @(t) 2+2*sin(2*pi*t);
LambdaReverse = @(t) 2+2*cos(2*pi*t);
L = @(t) 2+2*cos(2*pi*t+pi/6);
LReverse = @(t) 2+2*sin(2*pi*t+pi/6);
naturalBirdDeath = 0.3;
Cb = 3;
gamma = 0.1;
humansFearFlu = 5;
infectionRate = .1;
fluSurvivalRate = 0.9;

timesteps = 1000;

trackX = zeros(1,timesteps);
trackZ = zeros(1,timesteps);

trackYNsusc = zeros(1,timesteps);
trackYNinf = zeros(1,timesteps);
trackYDJsusc = zeros(1,timesteps);
trackYDJinf = zeros(1,timesteps);
trackYSsusc = zeros(1,timesteps);
trackYSinf = zeros(1,timesteps);

X = 1;
Z = 1;

YNsusc = 1;
YNinf = 0.1;
YDJsusc = 0.0;
YDJinf = 0.0;
YSsusc = 0.0;
YSinf = 0.0;

Xplus = X;
Zplus = Z;

YNsuscplus = YNsusc;
YNinfplus = YNinf;
YDJsuscplus = YDJsusc;
YDJinfplus = YDJinf;
YSsuscplus = YSsusc;
YSinfplus = YSinf;

for index=1:timesteps
    t = (index-1)*20/timesteps;

    trackX(index) = X;
    trackZ(index) = Z;

    trackYNsusc(index) = YNsusc;
    trackYNinf(index) = YNinf;
    trackYDJsusc(index) = YDJsusc;
    trackYDJinf(index) = YDJinf;
    trackYSsusc(index) = YSsusc;
    trackYSinf(index) = YSinf;

    if YNinf+YNsusc < 0.0001
        YNsusc = 0;
        YNinf = 0;
        NnewSusceptibleProportion = 1;
    else
        NnewSusceptibleProportion = exp(YNinf*(log(1-infectionRate))/(YNinf+YNsusc));
    end

    if YDJinf+YDJsusc < 0.0001
        YDJsusc = 0;
        YDJinf = 0;
        DJnewSusceptibleProportion = 1;
    else
        DJnewSusceptibleProportion = exp(YDJinf*(log(1-infectionRate))/(YDJinf+YDJsusc));
    end

    if YSinf+YSsusc < 0.0001
        YSsusc = 0;
        YSinf = 0;
        SnewSusceptibleProportion = 1;
    else
        SnewSusceptibleProportion = exp(YSinf*(log(1-infectionRate))/(YSinf+YSsusc));
    end

    % Susceptible birds north of Djoudj become infected
    YNsusc = YNsusc*NnewSusceptibleProportion;
    % Infected birds north of Djoudj die from the flu (within a month) 
    YNinf = YNsusc*(1-NnewSusceptibleProportion) + YNinf*fluSurvivalRate;
    YN = YNsusc + YNinf;

    % Susceptible Djoudj birds become infected
    YDJsusc = YDJsusc*DJnewSusceptibleProportion;
    % Susceptible Djoudj birds die from the flu (within a month)
    YDJinf = YDJsusc*(1-DJnewSusceptibleProportion) + YDJinf*fluSurvivalRate;
    YDJ = YDJsusc + YDJinf;

    % Susceptible birds south of Djoudj become infected
    YSsusc = YSsusc*SnewSusceptibleProportion;
    % Susceptible birds south of Djoudj die from the flu
    YSinf = YSsusc*(1-SnewSusceptibleProportion) + YSinf*fluSurvivalRate;
    YS = YSsusc + YSinf;

    XDeriv = r1*X*(1-X/k1) - alpha*X*YDJ/(X+D) - h1(t)*X*Z;

    % Susceptible birds north of Djoudj
    YNsuscDeriv = -naturalBirdDeath*YNsusc ... % die naturally
                    - Lambda(t)*YNsusc ... % migrate to Djoudj
                    + LambdaReverse(t)*YDJsusc ... % migrate from Djoudj
                    + r2*YNsusc*(1-YN/k2); % are born

    % Infected birds north of Djoudj
    YNinfDeriv = -naturalBirdDeath*YNinf ... % die naturally 
                    - Lambda(t)*YNinf ... % migrate to Djoudj (in 1-3 weeks)
                    + LambdaReverse(t)*YDJinf ... % migrate from Djoudj
                    + r2*YNinf*(1-YN/k2); % are born

    % Susceptible Djoudj birds
    YDJsuscDeriv = -naturalBirdDeath*YDJsusc ... % Die naturally
                    + beta*X*YDJsusc/(X+D) ... % Proliferate if there are fish
                    - h2(t)*YDJsusc ... % Are poached
                    - lambda(t)*YDJsusc*Z ... % Are hampered by tourists
                    + Lambda(t)*YNsusc ... % Migrate from the north
                    - L(t)*YDJsusc ... % Migrate south
                    + LReverse(t)*YSsusc ... % Migrate from the south
                    - LReverse(t)*YDJsusc ... % Migrate north
                    + r2*YDJsusc*(1-YDJ/k2); % Are born

    % Susceptible Djoudj birds
    YDJinfDeriv = -naturalBirdDeath*YDJinf ... % Die naturally
                    + beta*X*YDJinf/(X+D) ... % Proliferate if there are fish
                    - h2(t)*YDJinf ... % Are poached
                    - lambda(t)*YDJinf*Z ... % Are hampered by tourists
                    + Lambda(t)*YNinf ... % Migrate from the north
                    - L(t)*YDJinf ... % Migrate south
                    + LReverse(t)*YSinf ... % Migrate from the south
                    - LReverse(t)*YDJinf ... % Migrate north
                    + r2*YDJinf*(1-YDJ/k2); % Are born

    % Susceptible birds south of Djoudj
    YSsuscDeriv = -naturalBirdDeath*YSsusc ... % die naturally
                    + L(t)*YDJsusc ... % emigrate from Djoudj
                    - LReverse(t)*YSsusc ... % migrate to Djoudj
                    + r2*YSsusc*(1-YS/k2); % Are born

    % Susceptible birds south of Djoudj
    YSinfDeriv = -naturalBirdDeath*YSinf ... % die naturally
                    + L(t)*YDJinf ... % emigrate from Djoudj
                    - LReverse(t)*YSinf ...% migrate to Djoudj
                    + r2*YSinf*(1-YS/k2); % Are born

    % Tourism is high when there are healthy birds, low when there are sick
    % birds, and low when there's too many people
    ZDeriv = Cb*YDJsusc*Z/(Z+1) - gamma*Z^2 - humansFearFlu*YDJinf*Z;

    Xplus = X + XDeriv * 20/timesteps;
    YNsuscplus = YNsusc + YNsuscDeriv * 20/timesteps;
    YNinfplus = YNinf + YNinfDeriv * 20/timesteps;
    YDJsuscplus = YDJsusc + YDJsuscDeriv * 20/timesteps;
    YDJinfplus = YDJinf + YDJinfDeriv * 20/timesteps;
    YSsuscplus = YSsusc + YSsuscDeriv * 20/timesteps;
    YSinfplus = YSinf + YSinfDeriv * 20/timesteps;
    Zplus = Z + ZDeriv * 20/timesteps;

    if (Xplus < 0) X = 0; else X = Xplus; end
    if (YNsuscplus < 0) YNsusc = 0; else YNsusc = YNsuscplus; end
    if (YNinfplus < 0) YNinf = 0; else YNinf = YNinfplus; end
    if (YDJsuscplus < 0) YDJsusc = 0; else YDJsusc = YDJsuscplus; end
    if (YDJinfplus < 0) YDJinf = 0; else YDJinf = YDJinfplus; end
    if (YSsuscplus < 0) YSsusc = 0; else YSsusc = YSsuscplus; end
    if (YSinfplus < 0) YSinf = 0; else YSinf  = YSinfplus; end
    if (Zplus < 0) Z = 0; else Z = Zplus; end
end

% tiles = tiledlayout(5,1);
tiles = tiledlayout(3,1);

nexttile

t = linspace(0,20,timesteps); 
plot(t, trackX, '-');
title('Prey Population')
xlabel('Years');
ylabel('Population');
ylim([0,1]);

nexttile

plot(t, trackYNsusc+trackYNinf, '-');
hold on;
plot(t, trackYSsusc+trackYSinf, '-');
hold on;
plot(t, trackYDJsusc+trackYDJinf, '-');
title('Global Waterbird Population');
xlabel('Years');
ylabel('Population');
ylim([0,1]);
legend('Northern', 'Southern', 'Djoudj');


% plot(t, trackYNsusc+trackYNinf, '-');
% hold on;
% plot(t, trackYNsusc, '-');
% hold on;
% plot(t, trackYNinf);
% title('Northern Predator Population');
% xlabel('Years');
% ylabel('Population');
% ylim([0,5]);
% legend('Northern', 'Southern', 'Djoudj');

% nexttile
% 
% plot(t, trackYDJsusc+trackYDJinf, '-');
% hold on;
% plot(t, trackYDJsusc, '-');
% hold on;
% plot(t, trackYDJinf);
% title('Djoudj Predator Population');
% xlabel('Years');
% ylabel('Population');
% ylim([0,5]);
% 
% nexttile
% 
% plot(t, trackYSsusc+trackYSinf, '-');
% hold on;
% plot(t, trackYSsusc, '-');
% hold on;
% plot(t, trackYSinf);
% title('Southern Predator Population');
% xlabel('Years');
% ylabel('Population');
% ylim([0,5]);

nexttile

plot(t, trackZ, '- ');
title('Tourist Population');
xlabel('Years');
ylabel('Population');
ylim([0,2]);
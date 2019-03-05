% Read eikonal densities and elastic amplitudes
%
% mikael.mieskolainen@cern.ch, 2018

addpath ../../../../#matlabcodes
clear; close all;
figure;

while (true)

% Indices
Re_ind = 2;
Im_ind = 3;

sqrts  = [62 546 1800 7000 13000 60000];
s      = sqrts.^2;
factor = [20000 400 20 1 1/20 1/1000];

GeV2fm = 0.1973;
GeV2mb = 0.389;

% Differential cross section dsigma/dt
X = {};
legends = {};
for i = 1:length(sqrts)
    X{i} = csvread(sprintf('elastic_%0.0f.csv', sqrts(i)));
    if (sqrts(i) < 1000)
        legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i));
    elseif (sqrts(i) < 7000)
        legends{i} = sprintf('$\\sqrt{s} = %0.2f$ TeV', sqrts(i)/1e3);
    else
        legends{i} = sprintf('$\\sqrt{s} = %1.0f$ TeV', sqrts(i)/1e3);    
    end
end

% Generate dsigma/dt
dxs = {};
for i = 1:length(sqrts)
    dxs{i} = abs( X{i}(:,Re_ind) + 1i*X{i}(:,Im_ind) ).^2 / (16*pi*sqrts(i)^4) * GeV2mb;
end

% Font
for i = 1:length(sqrts)
    ind = find(X{i} == 3.0);
    if (factor(i) >= 1) 
        txt = sprintf('$\\times \\, %0.0f$', factor(i));
    elseif (factor(i) > 0.001)
        txt = sprintf('$\\times \\, %0.2f$', factor(i));    
    else
        txt = sprintf('$\\times \\, %0.3f$', factor(i));
    end
    text(X{i}(ind), factor(i)*dxs{i}(ind)*1.6, txt, 'interpreter','latex'); hold on;
end


for i = 1:length(sqrts)
    plot(X{i}(:,1), factor(i)*dxs{i}, '.', 'linewidth', 1.1); hold on;
end
set(gca,'yscale','log');
%set(gca,'xscale','log');
%axis tight;

% Read in TOTEM data
TOTEM_high_t = dlmread('../HEPdata/TOTEM_7.csv');
TOTEM_low_t  = dlmread('../HEPdata/TOTEM_7_low_t.csv');

TOTEM = [TOTEM_low_t; TOTEM_high_t];

errorbar(TOTEM(:,1), TOTEM(:,4) , ...
         sqrt(TOTEM(:,6).^2 + TOTEM(:,8).^2), ...
         sqrt(TOTEM(:,5).^2 + TOTEM(:,7).^2), 'k.', 'CapSize', 1);
%legends{length(legends)+1} = 'TOTEM $\sqrt{s} = 7$ TeV';

% Read in Fermilab data
ABE = dlmread('../HEPdata/ABE_1994_1800.csv');

errorbar(ABE(:,1), ABE(:,4)*factor(3), ...
         sqrt(ABE(:,6).^2 ), ...
         sqrt(ABE(:,5).^2 ), '.', 'color', [1 0 1]*0.3, 'CapSize', 1);
%legends{length(legends)+1} = 'D0 $\sqrt{s} = 1.96$ TeV';


% Read in FNAL data
FNAL = dlmread('../HEPdata/FNAL_546.csv');
SPS  = dlmread('../HEPdata/SPS_546.csv');
FNAL = [FNAL; SPS];

errorbar(FNAL(:,1), FNAL(:,4)*factor(2), ...
         sqrt(FNAL(:,6).^2 ), ...
         sqrt(FNAL(:,5).^2 ), '.', 'color', [1 1 0]*0.3, 'CapSize', 1);
%legends{length(legends)+1} = 'E741, SPS $\sqrt{s} = 546$ GeV';

%
% Read in ISR data
ISR = dlmread('../HEPdata/ISR_62.csv');

errorbar(ISR(:,1), ISR(:,4)*factor(1), ...
         sqrt(ISR(:,6).^2 ), ...
         sqrt(ISR(:,5).^2 ), '.', 'color', [1 1 0]*0.3, 'CapSize', 1);
%legends{length(legends)+1} = 'ISR $\sqrt{s} = 62.5$ GeV';
%}

% 3-gluon exchange

tval = linspace(1.5,8,1000);
scale = 3e-2; % Arbitrary scale
plot(tval, scale*tval.^(-8), 'k-.');

legends{length(legends)+1} = '$t^{-8}$';
%}

l = legend(legends); set(l,'interpreter','latex','fontsize', 8); legend('boxoff');
set(l,'interpreter','latex'); legend('boxoff');
xlabel('$-t$ (GeV$^2$)','interpreter','latex');
ylabel('$d\sigma/dt$ (mb/GeV$^2$)','interpreter','latex');
axis square;
axis([0 4.0 1e-10 1e6]);
xticks(linspace(0, 8, 17));

drawnow;
pause(3);
clf();
end

print -dpdf dsigmadt.pdf


%% Local B-slope

fig{1} = figure;
fig{2} = figure;

partial = {};

legends = {};
for i = 1:length(sqrts)
    legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i) );
end

for k = 1:2
    figure(fig{k});
    % Restart color order
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %plot(linspace(0,4,2), zeros(2,1), 'k-'); hold on;
    
    for i = 1:length(dxs)
        % d/dt ln(dsigma/dt)
        delta = X{i}(2,1) - X{i}(1,1);
        partial{i} = - diff(log(dxs{i})) / delta;
    
        plot(X{i}(1:end-1,1), partial{i}); hold on;
    end
    if (k == 1)
        axis([0 4 -inf inf]);
        set(gca,'xtick', linspace(0,4,11));
        l = legend(legends); set(l,'interpreter','latex'); legend('boxoff');
        set(l,'interpreter','latex'); legend('boxoff');    
    end
    if (k == 2)
        axis([0 0.25 10 35]);
        set(gca,'ytick', linspace(10,35,11));
    end
    
        xlabel('$-t$ (GeV$^2$)','interpreter','latex'); axis square;
        ylabel('$B \equiv \frac{d}{dt}\ln(d\sigma/dt)$','interpreter','latex');
end

%% B(t = 0)

figure;

B_t0 = zeros(length(sqrts),1);

for i = 1:length(sqrts)
   B_t0(i) = partial{i}(1);
end
plot(sqrts, B_t0, 'ks-');
set(gca,'xscale','log');

xlabel('$\sqrt{s}$ (GeV)','interpreter','latex');
ylabel('$B(t=0)$ (GeV$^{-2}$)','interpreter','latex');
axis square; axis tight;

%% Eikonal densities/opacities

figure;
X = {};
legends = {};
for i = 1:length(sqrts)
    X{i} = csvread(sprintf('density_%0.0f.csv', sqrts(i)));
    legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i));
end

for i = 1:length(X)
    plot(GeV2fm*X{i}(:,1), X{i}(:,Im_ind),'-'); hold on;
end

% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

for i = 1:length(X)
    plot(GeV2fm*X{i}(:,1), X{i}(:,Re_ind),'--'); hold on;
end
%set(gca,'yscale','log');
set(gca,'xscale','log');
axis tight;

l = legend(legends); set(l,'interpreter','latex'); legend('boxoff');
xlabel('$b$ (fm)','interpreter','latex');
ylabel('$\Omega(s,b)$','interpreter','latex'); axis square;


%% Elastic amplitudes in b-space

figure;
X = {};
legends = {};
for i = 1:length(sqrts)
    X{i} = csvread(sprintf('density_%0.0f.csv', sqrts(i)));
    legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i));
end
for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(X{i}(:,Re_ind) + 1i*X{i}(:,Im_ind))/2));
    plot(GeV2fm*X{i}(:,1), imag(A), '-'); hold on;
end
% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(X{i}(:,Re_ind) + 1i*X{i}(:,Im_ind))/2));
    plot(GeV2fm*X{i}(:,1), real(A), '--'); hold on;
end

%set(gca,'yscale','log');
%set(gca,'xscale','log');
%axis tight;
axis([0 4 0 1]); axis square;

l = legend(legends); set(l,'interpreter','latex'); legend('boxoff');
xlabel('$b$ (fm)','interpreter','latex');
ylabel('$A_{el}(s,b) = i(1 - \exp(i\Omega(s,b)/2))$','interpreter','latex');
set(gca,'xtick', linspace(0,4,11));


%% Real and Imag parts of elastic amplitudes in b-space

figure;
for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(X{i}(:,Re_ind) + 1i*X{i}(:,Im_ind))/2));
    plot(real(A), imag(A), '.'); hold on; axis square; axis tight;
    if i == length(sqrts)
        for k = 5:40:length(GeV2fm*X{i}(:,1))/9
            text(real(A(k)), imag(A(k)), sprintf('%0.2f fm', GeV2fm*X{i}(k,1)), 'interpreter','latex');
        end
    end
end

xlabel('Re [$A_{el}(s,b)$]','interpreter','latex');
ylabel('Im [$A_{el}(s,b)$]','interpreter','latex');


%% Real and Imag parts of elastic amplitudes in b-space

figure;
for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(X{i}(:,Re_ind) + 1i*X{i}(:,Im_ind))/2));
    plot(GeV2fm*X{i}(:,1), abs(real(A) ./ imag(A)), '-'); hold on; axis square; axis tight;
end
set(gca,'yscale','log');
xlabel('$b$ (fm)','interpreter','latex');
ylabel('$|$Re [$A_{el}(s,b)]/$ Im [$A_{el}(s,b)]|$','interpreter','latex');

axis([0 inf 1e-2 10]);


%% Elastic amplitudes in t-space

figure;
X = {};
legends = {};
for i = 1:length(sqrts)
    X{i} = csvread(sprintf('elastic_%0.0f.csv', sqrts(i)));
    legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i));
end

for i = 1:length(sqrts)
    plot(X{i}(:,1), X{i}(:,Im_ind) / s(i), '-'); hold on;
end
% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;
for i = 1:length(sqrts)
    plot(X{i}(:,1), X{i}(:,Re_ind) / s(i), '--'); hold on;
end

%set(gca,'yscale','log');
set(gca,'xscale','log');
axis tight; axis square;
axis([0 inf -20 inf]);

l = legend(legends); set(l,'interpreter','latex'); legend('boxoff');
xlabel('$-t$ (GeV$^2$)','interpreter','latex');
ylabel('$A_{el}(s,t) / s$','interpreter','latex');



%% Real and Imag parts of elastic amplitudes in t-space

figure;
for i = 1:length(sqrts)
    A = (X{i}(:,Re_ind) + 1i*X{i}(:,Im_ind))/ s(i);
    plot(real(A), imag(A), '-'); hold on; axis square; axis tight;
    if i == length(sqrts)
        for k = 1:45:length(GeV2fm*X{i}(:,1))/8
            text(real(A(k)), imag(A(k)), sprintf('%0.2f', -X{i}(k,1)), 'interpreter','latex');
        end
    end
end
axis([-5 5 -10 100]);
xlabel('Re [$A_{el}(s,t) / s$]','interpreter','latex');
ylabel('Im [$A_{el}(s,t) / s$]','interpreter','latex');
set(gca,'xtick',linspace(-5,5,11));


%% Re{A_el(s,t)}/Im{A_el(s,t)}

figure;
for i = 1:length(sqrts)
    A = (X{i}(:,Re_ind) + 1i*X{i}(:,Im_ind))/ s(i);
    plot(X{i}(:,1), abs(real(A) ./ imag(A)), '-'); hold on;
end
axis square; axis tight;
set(gca,'yscale','log');

xlabel('$-t$ (GeV$^2$)','interpreter','latex');
ylabel('$|$Re [$A_{el}(s,t)]/$ Im [$A_{el}(s,t)]|$','interpreter','latex');

axis([0 inf 1e-3 20]);


%% Rho-parameter rho = Re{A_el(s,t=0)}/Im{A_el(s,t=0)}

figure;
X = {};
for i = 1:length(sqrts)
    X{i} = csvread(sprintf('elastic_%0.0f.csv', sqrts(i)));
end

% Construct rho
rho = zeros(length(X),1);
for i = 1:length(sqrts)
    rho(i) = X{i}(1,Re_ind) / X{i}(1,Im_ind);
end
plot(sqrts, rho, 'ks-'); hold on;

%set(gca,'yscale','log');
set(gca,'xscale','log');
axis tight; axis square;

xlabel('$\sqrt{s}$ (GeV)','interpreter','latex');
ylabel('$\rho \equiv$ Re [$A_{el}(s,t=0)]/$ Im [$A_{el}(s,t=0)$]','interpreter','latex');

axis([0 inf 0.09 0.14]);



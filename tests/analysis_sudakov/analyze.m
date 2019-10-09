% Read out Sudakov supression integral arrays
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath ../mcodes

sqrts = 13000;
pdf = {'CT10nlo'};
%pdf = {'MSTW2008lo68cl'};
%pdf = {'MMHT2014lo68cl'};

legends = {};
for i = 1:length(pdf)
    % Read the Sudakov array file
    [~,filename] = system(sprintf('ls ../../sudakov/SUDA_%0.0f_%s_*', sqrts, pdf{i}));
    filename = filename(1:end-1); % Remove \n
    T{i} = csvread(filename);
    
    % Read the Shuvaev array file
    [~,filename] = system(sprintf('ls ../../sudakov/SHUV_%0.0f_%s_*', sqrts, pdf{i}));
    filename = filename(1:end-1); % Remove \n
    H{i} = csvread(filename);
    
    %legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i));
end

PRINT_ON = true;

%% Shuvaev pdf array H(q2, log(x))

X = H{1};

close all;

% These need to match the array file!
Nq2  = 800;
Nlnx = 800;

k = 1;
Hmatrix = zeros(Nq2, Nlnx);
q2val   = zeros(Nq2,1);

for z = 1:Nq2
    start    = (z-1)*(Nlnx+1) + 1;
    stop     = start + Nlnx - 1;
    
    Hval = X(start:stop, 3);
    
    Hmatrix(k,:)   = Hval;
    q2val(k) = X(start,1);
    
    k = k + 1;
end
lnxval = X(1:Nlnx,2);
Hmatrix = Hmatrix';

imagesc(q2val, exp(lnxval), Hmatrix);
set(gca,'YDir','Normal');
xlabel('$Q^2$ (GeV$^2$)','interpreter','latex');
ylabel('$x$',    'interpreter','latex');
title('$H(Q^2,x)$',    'interpreter','latex');
colorbar;
%caxis([0 1]);
axis([0 max(q2val) 0 1]);
axis square;

%colormap bone
shading interp


%% Sudakov array T(q2, log(M))

X = T{1};

close all;

% These need to match the array file!
Nq2  = 800;
NlnM = 800;

legs = {};
k = 1;
Tmatrix = zeros(Nq2, NlnM);
q2val   = zeros(Nq2,1);

for z = 1:Nq2
    start    = (z-1)*(NlnM+1) + 1;
    stop     = start + NlnM - 1;
    
    Tval = X(start:stop, 3);
    
    Tmatrix(k,:)   = Tval;
    q2val(k) = X(start,1);
    
    k = k + 1;
end
lnMval = X(1:NlnM,2);

Tmatrix = Tmatrix';
size(Tmatrix)

% 2D
figure;
imagesc(q2val, exp(lnMval), log10(Tmatrix)); hold on;

%[X,Y] = meshgrid(q2val, exp(lnMval));
%contour(X,Y,log10(Tmatrix + 1), 5);

set(gca,'YDir','Normal');
xlabel('$Q_t^2$ (GeV$^2$)','interpreter','latex');
ylabel('$\mu$ (GeV)',    'interpreter','latex');
%title('$T_g(Q_t^2,\mu^2)$',    'interpreter','latex');
colorbar;
h = colorbar;
ylabel(h, '$T(Q_t^2,\mu^2)$','interpreter','latex');

caxis([-3 0]);
axis square;

colormap bone(1000)
shading interp

% PRINT OUT
if (PRINT_ON)
filename = sprintf('./figs/sudakov_2D.pdf');
print(gcf, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));
end

%% Sudakov 1D over mu
fig2 = figure;
legs = {};
steps = [1 3 10 31 96 295 800];

for k = steps
plot(exp(lnMval), Tmatrix(:,k)); hold on;
set(gca,'yscale','log');
set(gca,'xscale','log');

Q2 = q2val(k);
if (Q2 < 100)
legs{end+1} = sprintf('$Q_t^2 = %0.2g$ GeV$^2$', Q2);
else
legs{end+1} = sprintf('$Q_t^2 = %0.0f$ GeV$^2$', Q2);    
end
end
axis([0 1200 0.99e-4 1.5]); axis square;

l = legend(legs); set(l,'interpreter','latex','location','southeast');
xlabel('$\mu$ (GeV)', 'interpreter','latex');
ylabel('$T(Q_t^2, \mu^2)$', 'interpreter','latex');

% PRINT OUT
if (PRINT_ON)
filename = sprintf('./figs/sudakov_1D_mu.pdf');
print(fig2, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''2 2 2 2'' %s %s', filename, filename));
end

%% Sudakov 1D over Q^2
close all;

fig3 = figure;
legs = {};
steps = [1, 75, 168, 250, 334, 410, 501];

for k = steps
plot(q2val, Tmatrix(k,:)); hold on;
set(gca,'yscale','log');
set(gca,'xscale','log');

mu = exp(lnMval(k));
if (mu < 100)
legs{end+1} = sprintf('$\\mu = %0.2g$ GeV', mu);
else
legs{end+1} = sprintf('$\\mu = %0.0f$ GeV', mu);
end

end
axis([0 max(q2val) 0.99e-4 1.5]);
axis square;

l = legend(legs); set(l,'interpreter','latex','location','southeast');
xlabel('$Q_t^2$ (GeV$^2$)', 'interpreter','latex');
ylabel('$T(Q_t^2, \mu^2)$', 'interpreter','latex');

% PRINT OUT
if (PRINT_ON)
filename = sprintf('./figs/sudakov_1D_Q2.pdf');
print(fig3, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''2 2 2 2'' %s %s', filename, filename));
end


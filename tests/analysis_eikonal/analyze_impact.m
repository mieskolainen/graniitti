% Read eikonal densities in impact parameter space
%
% mikael.mieskolainen@cern.ch, 2019

close all; clear;

addpath ../mcodes

% Run this first
setup;

fig0 = figure;

Xb = {};
legends = {};
for i = 1:length(sqrts)
    % Read the density array file
    [~,filename] = system(sprintf('ls ../../eikonal/MBT_2212_2212_%0.0f_*', sqrts(i)));
    filename = filename(1:end-1); % Remove \n
    Xb{i} = csvread(filename);
    
    legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i));
end


%% b-space profile (OMEGA)

fig0 = figure;

% undersampling for visualization (otherwise no-vector graphics)
step = 2;

% IMAG
for i = 1:length(Xb)
    plot(GeV2fm*Xb{i}(1:step:end,1), Xb{i}(1:step:end,Im_ind),'-'); hold on;
end

% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

% REAL
for i = 1:length(Xb)
    plot(GeV2fm*Xb{i}(1:step:end,1), Xb{i}(1:step:end,Re_ind),'--'); hold on;
end
%set(gca,'yscale','log');
%set(gca,'xscale','log');
axis tight;
axis([1e-2 2.5 0 8]); axis square;

title('Re[$\Omega$] (dashed), Im[$\Omega$] (solid)','interpreter','latex');

l = legend(legends); set(l,'interpreter','latex'); legend('boxoff');
xlabel('$b$ (fm)','interpreter','latex');
ylabel('$\Omega(s,b)$','interpreter','latex'); axis square;

% PRINT OUT
filename = sprintf('./figs/b_space_omega.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));

%% b-space profile ratio (OMEGA)

fig0 = figure;

% undersampling for visualization (otherwise no-vector graphics)
step = 2;

% ratio
for i = 1:length(Xb)
    plot(GeV2fm*Xb{i}(1:step:end,1), Xb{i}(1:step:end,Re_ind) ./ Xb{i}(1:step:end,Im_ind),'-'); hold on;
end

%set(gca,'yscale','log');
%set(gca,'xscale','log');
axis tight;
axis([0 4 0.12 0.24]); axis square;

l = legend(legends); set(l,'interpreter','latex','location','southeast'); legend('boxoff');
xlabel('$b$ (fm)','interpreter','latex');
ylabel('Re $[\Omega(s,b)]$ / Im $[\Omega(s,b)]$','interpreter','latex'); axis square;

% PRINT OUT
filename = sprintf('./figs/b_space_omega_ratio.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


%% Elastic amplitudes in b-space

fig0 = figure;

Xb = {};
legends = {};
for i = 1:length(sqrts)
    
    % Read the density array file
    [~,filename] = system(sprintf('ls ../../eikonal/MBT_2212_2212_%0.0f_*', sqrts(i)));
    filename = filename(1:end-1); % Remove \n
    Xb{i} = csvread(filename);
    
    legends{i} = sprintf('$\\sqrt{s} = %0.0f$ GeV', sqrts(i));
end

% Undersampling for graphics reasons
step = 2;
for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(Xb{i}(1:step:end,Re_ind) + 1i*Xb{i}(1:step:end,Im_ind))/2));
    plot(GeV2fm*Xb{i}(1:step:end,1), imag(A), '-'); hold on;
end
% Restart color order
ax = gca;
ax.ColorOrderIndex = 1;

for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(Xb{i}(:,Re_ind) + 1i*Xb{i}(:,Im_ind))/2));
    plot(GeV2fm*Xb{i}(:,1), real(A), '--'); hold on;
end

title('Re[A] (dashed), Im[A] (solid)','interpreter','latex');
%set(gca,'yscale','log');
%set(gca,'xscale','log');
%axis tight;
axis([0 4 0 1]); axis square;

l = legend(legends); set(l,'interpreter','latex'); legend('boxoff');
xlabel('$b$ (fm)','interpreter','latex');
ylabel('$A_{el}(s,b) = i[1 - \exp(i\Omega(s,b)/2)]$','interpreter','latex');
set(gca,'xtick', linspace(0,4,11));

% PRINT OUT
filename = sprintf('./figs/b_space_amplitude.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


%% Real and Imag parts of elastic amplitudes in b-space (spiral plot)

fig0 = figure;
for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(Xb{i}(:,Re_ind) + 1i*Xb{i}(:,Im_ind))/2));
    plot(real(A), imag(A), '.'); hold on; axis square; axis tight;
    if i == length(sqrts)
        for k = 5:40:length(GeV2fm*Xb{i}(:,1))/9
            text(real(A(k)), imag(A(k)), sprintf('%0.2f fm', GeV2fm*Xb{i}(k,1)), 'interpreter','latex');
        end
    end
end

xlabel('Re [$A_{el}(s,b)$]','interpreter','latex');
ylabel('Im [$A_{el}(s,b)$]','interpreter','latex');



%% Real and Imag parts of elastic amplitudes in b-space

fig0 = figure;
for i = 1:length(sqrts)
    A = 1i*(1 - exp(1i*(Xb{i}(:,Re_ind) + 1i*Xb{i}(:,Im_ind))/2));
    plot(GeV2fm*Xb{i}(:,1), real(A) ./ imag(A), '-'); hold on; axis square; axis tight;
end
set(gca,'yscale','log');
xlabel('$b$ (fm)','interpreter','latex');
ylabel('Re [$A_{el}(s,b)]$ / Im [$A_{el}(s,b)]$','interpreter','latex');

axis([0 inf 1e-2 1]);

l = legend(legends); set(l,'interpreter','latex'); legend('boxoff');
legend('location','northwest');


% PRINT OUT
filename = sprintf('./figs/b_amplitude_re_im.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));

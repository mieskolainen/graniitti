% Cross sections as a function of energy, screening on/off
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath /home/user/cernbox/#matlabcodes

% Copy this to the bash script looping over CMS energies
sqrts_values = logspace(1.3, 6.3, 10);
fprintf('E=');
for i = 1:length(sqrts_values)
    if (i < length(sqrts_values))
        marker = ',';
    else
        marker = '';
    end
    fprintf('%0.6E%s', sqrts_values(i), marker);
end
fprintf('\n\n');


%% SET MC and DATA here

% 1=sqrts 2=total 3=inel 4=el 5=xs0 6=xs1 7=xs2  [barn]
MC_S_false = dlmread('scan_screening_false.csv', '\t', 1, 0);
MC_S_true  = dlmread('scan_screening_true.csv',  '\t', 1, 0);

% Take non-screened
MC = MC_S_false;

% Take screened
MC = [MC MC_S_true(:,5:7)];

% Scale x-sections to ubarn
MC(:,2:end) = MC(:,2:end) * 1e6;


%% Plot cross sections
close all;

colors = {'k-','k--','k:', 'r-','r--','r:'};
kk = 1;
h = {};
for i = 5:10
    h{end+1} = plot(MC(:,1), MC(:,i), colors{kk}); hold on;
    kk = kk + 1;
end

l = legend('$\pi^+\pi^-_{EL}$','$\pi^+\pi^-_{SD}$','$\pi^+\pi^-_{DD}$', ...
           '$\pi^+\pi^-_{EL}$ ($S^2$)','$\pi^+\pi^-_{SD}$ ($S^2$)','$\pi^+\pi^-_{DD}$ ($S^2$)');
set(l,'interpreter','latex','location','northwest'); legend('boxoff');

set(gca,'xscale','log');
set(gca,'yscale','log');

axis tight;

xlabel('$\sqrt{s}$ (GeV)','interpreter','latex');
ylabel('$\sigma$ ($\mu$b)','interpreter','latex');

axis([10 max(MC(:,1)) 1e-1 1e2]);

% PRINT OUT
filename = sprintf('./figs/xsec.pdf');
print(gcf, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


%% Plot ratios
close all;

colors = {'k-','k--','k:', 'r-','r--','r:'};
kk = 1;
h = {};

pairs = [5 8;
         6 9;
         7 10];

for i = 1:3
    h{end+1} = plot(MC(:,1), MC(:,pairs(i,2)) ./ MC(:,pairs(i,1)), colors{kk}); hold on;
    kk = kk + 1;
end

l = legend('$\pi^+\pi^-_{EL}$','$\pi^+\pi^-_{SD}$','$\pi^+\pi^-_{DD}$');

set(l,'interpreter','latex','location','northeast'); legend('boxoff');

set(gca,'xscale','log');
%set(gca,'yscale','log');

axis tight;

xlabel('$\sqrt{s}$ (GeV)','interpreter','latex');
ylabel('$\langle S^2 \rangle$','interpreter','latex');

axis([10 max(MC(:,1)) 0 0.5]);

% PRINT OUT
filename = sprintf('./figs/xsecratios.pdf');
print(gcf, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


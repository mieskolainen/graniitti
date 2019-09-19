% Diffractive cross sections as a function of energy
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

addpath /home/user/cernbox/#matlabcodes

% Copy this to the bash script looping over CMS energies
sqrts_values = logspace(1.3, 6.3, 100);
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

% Indices
tot_id = 2;
inel_id = 3;
el_id = 4;
sd_id = 7; % 7
dd_id = 8; % 8

% IMPORTANT!!! Set here MANUALLY the
% Triple Pomeron Coupling used in the calculation -> otherwise the re-fit
% values make no sense!
g3P_orig = 0.18;

% Coupling grid scan values
P3values = linspace(0.05, 0.25, 1e3);

% sqrts total inel el xs0(sd) xs1(dd) xs2(sd & xi < 0.05) xs3(dd & xi < 0.05)  [barn]
MC = dlmread('scan.csv', '\t', 1, 0);

% Scale x-sections to mb
MC(:,2:end) = MC(:,2:end) * 1e3;

colors = {'k-','r-','k--'};
kk = 1;

h = {};
for i = [tot_id inel_id el_id]
    h{end+1} = plot(MC(:,1), MC(:,i), colors{kk}); hold on;
    kk = kk + 1;
end

% ************************************************************************
% TOTAL cross section measurements

% TOTEM Run 1
X{1} = [7000 98.3 sqrt(0.2^2 + 2.8^2);             % Total
        7000 73.5 sqrt(0.6^2 + ((1.8--1.3)/2)^2);  % Inelastic
        7000 24.8 sqrt(0.2^2 + 1.2^2)];            % Elastic

% CDF
X{2} = [1800 80.03 2.24; % Total
        1800 60.6 sqrt(2.24^2 + 0.85^2); % INELASTIC (SELF CALCULATED FROM TOTAL & EL)
        1800 19.7 0.85];  % Elastic

% ISR
X{3} = [62 43.9 0.3;   % Total
        62 36.29 sqrt(0.3^2 + 0.22^2); % Self calculated inelastic
        62 7.61 0.22]; % Elastic
    
X{4} = [53 42.9 0.2;   % Total
        53 35.53 sqrt(0.2^2 + 0.13^2); % Self calculated inelastic
        53 7.37 0.13]; % Elastic
    
X{5} = [31 40.2 0.2;   % Total
        31 33.08 sqrt(0.2^2 + 0.34^2); % Self calculated inelastic
        31 7.12 0.34]; % Elastic

marks = {'k.','ko','ks','ks','ks'};
h = {};
for i = 1:length(X)
    h{end+1} = plot(X{i}(:,1), X{i}(:,2), marks{i});
    errorbar(X{i}(:,1), X{i}(:,2), X{i}(:,3), marks{i}); hold on;
end

l = legend([h{[1:3]}], ...
{'LHC (TOTEM)','Tevatron (CDF)','ISR'} );
set(l,'interpreter','latex','location','northwest'); legend('boxoff');

set(gca,'xscale','log');
%set(gca,'yscale','log');

%l = legend('$\sigma_{tot}$','$\sigma_{in}$','$\sigma_{el}$','$\sigma_{sd}$','$\sigma_{dd}$');
%set(l,'interpreter','latex','location','northwest'); legend('boxoff');
axis tight;
%axis square;

xlabel('$\sqrt{s}$ (GeV)','interpreter','latex');
ylabel('$\sigma$ (mb)','interpreter','latex');

axis([10 max(MC(:,1)) 0 110]);

% PRINT OUT
filename = sprintf('./figs/total.pdf');
print(gcf, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


%% SD Measurements

% [sqrts(s) value error]

XSD = {};

% -- ISR --
% J.C.M. Armitage et al., Nucl. Phys. B194, 365 (1982)

XSD{1} = [23.45 6.5 0.2;
        62.29 7.5 0.3];

% -- SppS --
% D.  Bernard  et  al.  (UA4  Collaboration),  Phys.  Lett.  B186,  227 (1987)
XSD{2} = [546 9.4 0.7];

% -- Tevatron --

% N.A. Amos  et al.  (E710  Collaboration), Phys.  Lett.  B 301 ,  313 (1993)
XSD{3} = [1800 9.4 1.4];

% CDF, https://journals.aps.org/prd/pdf/10.1103/PhysRevD.50.5535
XSD{4} = [546 7.89 0.33;
        1800 9.46 0.44];

% -- LHC --

% ALICE Run1
XSD{5} = [900  11.2  (1.6 - -2.1)/2
        2760 12.12 (3.9 - -5.3)/2
        7000 14.9  (3.4 - -5.9)/2];
    
% CMS Run1, https://arxiv.org/pdf/1503.08689.pdf
XSD{6} = [7000 8.84 sqrt(0.08^2 + 1.49^2 + 1.17^2)];


%% DD measurements

% [sqrts(s) value error]

XDD = {};

% -- ISR --
% 

% -- SppS --
% UA5
XDD{1} = [200 4.9 3];

% -- Tevatron --
% CDF
XDD{2} = [546  4 1.8;
        1800 6 1.5];

% -- LHC --
% ALICE Run1
XDD{3} = [900  5.6  2.0;
        2760 7.8  3.2;
        7000 9.0  2.6];

XDD{4} = [7000 5.17 sqrt(0.08^2 + 0.57^2 + 1.62^2)];    

%% 3P coupling fit

fig0 = figure;
chi2_values = zeros(length(P3values),1);

for k = 1:length(P3values)    
    
    % CHI^2
    chi2 = 0;
%
    % SD
    for i = 1:length(XSD)
        M = XSD{i};
        % Loop over measurements
        for j = 1:size(M,1)
            
            energy = M(j,1);
            % Find the closest energy value of MC
            [value,ind] = min(abs(energy - MC(:,1)));
            
            % NEW MC value
            MC_value = MC(ind, sd_id) / g3P_orig * P3values(k);
            
            cost =  ((MC_value - M(j,2))/M(j,3))^2;
            chi2 = chi2 + cost;
        end
    end
%}
    % DD
    for i = 1:length(XDD)
        M = XDD{i};
        % Loop over measurements
        for j = 1:size(M,1)
            
            energy = M(j,1);
            
            % Find the closest energy value of MC
            [value,ind] = min(abs(energy - MC(:,1)));
            
            % NEW MC value
            MC_value = (MC(ind, dd_id) / (g3P_orig^2)) * P3values(k)^2;
            
            cost =  ((MC_value - M(j,2))/M(j,3))^2;
            chi2 = chi2 + cost;
        end
    end
    chi2_values(k) = chi2;
    
    fprintf('chi2 = %0.5f, g3P = %0.5f \n', chi2, P3values(k));
end

% Find chi^2 + 1 (approx correspond to 1 sigma errors)
[minchi2,ind] = min(chi2_values);
lowind = 1;
for i = ind:-1:1
    if (chi2_values(i) > minchi2 + 1)
        lowind = i;
        break; 
    end
end
upind = 0;
for i = ind:1:length(chi2_values)
    if (chi2_values(i) > minchi2 + 1)
        upind = i;
        break; 
    end
end

g3P      = P3values(ind);
g3P_up   = P3values(upind);
g3P_down = P3values(lowind);

plot( ones(10,1)*P3values(lowind), linspace(0.5, max(chi2_values), 10), 'r--'); hold on;
plot( ones(10,1)*P3values(upind),  linspace(0.5, max(chi2_values), 10), 'r--');
plot(P3values, chi2_values);

xlabel('$g_{3P}$','interpreter','latex');
ylabel('$\chi^2$','interpreter','latex');
set(gca,'yscale','log');

axis([min(P3values) max(P3values) 0.5 100]);

fprintf('g3P = %0.3f (%0.3f ... %0.3f) \n', P3values(ind), P3values(lowind), P3values(upind));

l = legend('$\chi^2_{min} + 1$');
set(l,'interpreter','latex');

% PRINT OUT
filename = sprintf('./figs/g3p_chi2fit.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


%% RATIOPLOT

% Triple Pomeron Coupling variation (linear effect on SD cross section)
SD_up  = (MC(:,sd_id) / g3P_orig) * g3P_up;
SD     = (MC(:,sd_id) / g3P_orig) * g3P;
SD_low = (MC(:,sd_id) / g3P_orig) * g3P_down;

% Triple Pomeron Coupling variation (squared on DD cross section)
DD_up  = (MC(:,dd_id) / g3P_orig^2) * g3P_up^2;
DD     = (MC(:,dd_id) / g3P_orig^2) * g3P^2;
DD_low = (MC(:,dd_id) / g3P_orig^2) * g3P_down^2;

% Ratioplots

fig0 = figure;

% DIFF/TOT
plotfill(MC(:,1), sum(MC(:,el_id) + SD_up + DD_up, 2) ./ MC(:,tot_id), ...
         sum(MC(:,el_id) + SD_low + DD_low, 2) ./ MC(:,tot_id), [1 1 1]*0.7, [1 1 1]*0.5, 0.1); hold on;

h{1} = plot(MC(:,1), sum(MC(:,el_id) + SD + DD,2) ./ MC(:,tot_id), 'k-','linewidth',1.1); hold on;

% SD/INEL
plotfill(MC(:,1), SD_up ./ MC(:,inel_id), SD_low ./ MC(:,inel_id), [1 1 1]*0.7, [1 1 1]*0.5, 0.1); hold on;
h{2} = plot(MC(:,1), SD ./ MC(:,inel_id), 'r--');

% DD/INEL
plotfill(MC(:,1), DD_up ./ MC(:,inel_id), DD_low ./ MC(:,inel_id), [1 1 1]*0.7, [1 1 1]*0.5, 0.1); hold on;
h{3} = plot(MC(:,1), DD ./ MC(:,inel_id), 'k-.');

% EL/INEL
h{4} = plot(MC(:,1), MC(:,el_id) ./ MC(:,inel_id), 'k:');

% Pumplin bound
plot(MC(:,1), ones(length(MC(:,1)), 1)*0.5, 'r-','linewidth',1.1);
set(gca,'xscale','log');
set(gca,'yscale','log');

l = legend([h{1:4}], '$\sigma_{EL+SD+DD}/\sigma_{TOT}$', ...
    '$\sigma_{SD}/\sigma_{INEL}$', '$\sigma_{DD}/\sigma_{INEL}$', '$\sigma_{EL}/\sigma_{INEL}$', 'MP bound');
set(l,'interpreter','latex','location','southwest'); legend('boxoff');

xlabel('$\sqrt{s}$ (GeV)','interpreter','latex');
ylabel('','interpreter','latex');

axis([10 max(MC(:,1)) 1e-2 1]);

% PRINT OUT
filename = sprintf('./figs/ratio.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


%% SD PLOT

fig0 = figure;

plotfill(MC(:,1), SD_up, SD_low, [1 1 1]*0.7, [1 1 1]*0.5, 0.1); hold on;
h0 = plot(MC(:,1), SD, 'k-');

marks = {'k.','k^','ks','ko','rs','kx'};
h = {};
for i = 1:length(XSD)
    h{i} = plot(XSD{i}(:,1), XSD{i}(:,2), marks{i});
    errorbar(XSD{i}(:,1), XSD{i}(:,2), XSD{i}(:,3), marks{i}); hold on;
end

l = legend([h{1:length(XSD)}], ...
{'ISR $(M_X^2 < 0.05s)$','UA4, $M_X^2 < 0.05s$ (Sp$\bar{p}$S)', ...
'E710, $2$ GeV$^2 < M_X^2 < 0.05s$','CDF','ALICE, $M_X < 200$ GeV','CMS, $M_X^2 < 0.05s$'} );
set(l,'interpreter','latex','location','northwest'); legend('boxoff');

xlabel('$\sqrt{s}$ (GeV)',  'interpreter','latex');
ylabel('$\sigma_{SD}$ (mb)','interpreter','latex');

set(gca,'xscale','log');
axis([10 max(MC(:,1)) 0 22]);

% PRINT OUT
filename = sprintf('./figs/sd.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));


%% DD PLOT

fig0 = figure;

h0 = ...
plotfill(MC(:,1), DD_up, DD_low, [1 1 1]*0.7, [1 1 1]*0.5, 0.1); hold on;
h0 = plot(MC(:,1), DD, 'k-');

marks = {'ks','ko','rs','kx'};
h = {};
for i = 1:length(XDD)
    h{i} = plot(XDD{i}(:,1), XDD{i}(:,2), marks{i});
    errorbar(XDD{i}(:,1), XDD{i}(:,2), XDD{i}(:,3), marks{i}); hold on;
end

l = legend([h{1:length(XDD)}], {'UA5', 'CDF', 'ALICE', 'CMS'} );
set(l,'interpreter','latex','location','northwest'); legend('boxoff');

xlabel('$\sqrt{s}$ (GeV)','interpreter','latex');
ylabel('$\sigma_{DD}$ (mb)','interpreter','latex');

set(gca,'xscale','log');
axis([10 max(MC(:,1)) 0 14]);

% PRINT OUT
filename = sprintf('./figs/dd.pdf');
print(fig0, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''10 10 10 10'' %s %s', filename, filename));

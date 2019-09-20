% Read out Sudakov supression integral arrays
%
% mikael.mieskolainen@cern.ch, 2019
clear; close all;

sqrts = 13000;
pdf = {'MSTW2008lo68cl'};
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


%% Shuvaev pdf array H(q2, log(x))

X = H{1};

close all;

Nq2  = 300;
Nlnx = 3000;

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
title('log$_{10} \, H(Q^2, x)$',    'interpreter','latex');
colorbar;
%caxis([0 1]);
axis([0 max(q2val) 0 1]);
axis square;

%colormap bone
shading interp


%% Sudakov array T(q2, log(M))

X = T{1};

close all;

Nq2  = 300;
NlnM = 3000;

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

imagesc(q2val, exp(lnMval), Tmatrix);
set(gca,'YDir','Normal');
xlabel('$Q^2$ (GeV$^2$)','interpreter','latex');
ylabel('$\mu$ (GeV)',    'interpreter','latex');
title('log$_{10} \, T(Q^2,\mu)$',    'interpreter','latex');
colorbar;
caxis([0 1]);
axis square;

colormap summer(100)
shading interp


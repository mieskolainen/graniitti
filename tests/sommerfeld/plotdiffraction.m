% Plot Sommerfeld scalar diffraction simulation
%
% mikael.mieskolainen@cern.ch, 2019

clear;

A = csvread('../../3D.ascii');
X = zeros(max(A(:,1))+1, max(A(:,2))+1, max(A(:,3))+1, max(A(:,4))+1);

for k = 1:size(A,1)
    % +1 for matlab indexing
    ind = sub2ind(size(X), A(k,1) + 1, A(k,2) + 1, A(k,3) + 1, A(k,4) + 1);
    X(ind) = A(k,end-1) + 1i*A(k,end); 
end

%%
close all;

xval = linspace(-1,1,size(X,1));
yval = linspace(-1,1,size(X,2));
zval = linspace(-1,5.0,size(X,3));
tval = linspace( 0,1,size(X,4));

bound = max(abs(min(X(:))), abs(max(X(:))));
zlimits = [-bound, bound];

% t-list
tlist = 1:size(X,4);
tlist = 1;

for rep = 1:1 % Infinite loop
for t = tlist
    
    M = real( squeeze(X(1,:,:,t)) );
    imagesc(zval, yval, M);
    xlabel('$z$','interpreter','latex');
    ylabel('$x$','interpreter','latex');
    
    caxis(zlimits); %colorbar;
    hold on;
    
    plot([0 0], [-1.0 -0.8], '-k', 'linewidth', 2.0);
    plot([0 0], [-0.6 -0.2], '-k', 'linewidth', 2.0);
    plot([0 0], [ 0.2  0.6], '-k', 'linewidth', 2.0);
    plot([0 0], [ 0.8  1.0], '-k', 'linewidth', 2.0);
    
    set(gca,'YDir','normal');
    xticks([min(zval):0.5:max(zval)]);
    yticks([-1:0.2:1]);
    
    axis square;
    colormap hot;
    drawnow;
    pause(0.01);
end
end

% PRINT OUT
filename = sprintf('./figs/sommerfeld.pdf');
print(gcf, '-dpdf', filename);
system(sprintf('pdfcrop --margins ''2 2 2 2'' %s %s', filename, filename));

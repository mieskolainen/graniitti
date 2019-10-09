% Double ratioplot visualization
% 
% Output:  fig = Figure handle
%           ax = Two axes handles
%
% Usage example:
% 
%   [fig,ax] = ratioplot();
%   
%   axes(ax{1});
%   plot(x,y);
%   set(ax{1}, 'YScale', 'log');
%   axes(ax{2});
%   plot(x,y1./y2);
%
% mikael.mieskolainen@cern.ch, 2019

function [fig, ax] = ratioplot()

% Plot the results
fig = figure('Visible', 'Off', 'Units', 'Normalized', 'OuterPosition', [0 0 0.6 0.6]);

% Create upper ax
posvec1 = [0.13 0.273 0.775 0.662];
ax{1} = axes('Parent', fig, 'Position', posvec1);
hold(ax{1}, 'on');

% Box on
box(ax{1}, 'on');

% Set the remaining ax properties
set(ax{1},  'XMinorTick', 'on', 'XTick', [], 'YMinorTick', 'on');
hold(ax{1}, 'on');

% PLOT1 : Main plot
axis square; 

% Remove lowest tick mark, because it can overlap with the lower plot
%yTick = get(ax{1}, 'YTick');
%set(ax{1}, 'YTick', yTick(2:end));

axis normal;

% ------------------------------------------------------------------------
% PLOT2 : Ratio plot
% Create ax
screensize = get(0, 'screensize');
if (screensize(3) == 1920 && screensize(4) == 1200)
    posvec2 = [0.3372 0.132 0.36111 0.136];
else
    posvec2 = [0.362 0.132 0.30965 0.136]; % lenovo
end
ax{2} = axes('Parent', fig, 'Position', posvec2);

% Now new position for upper (need for exact horizontal alingment)
set(ax{1},'position', [posvec2(1) posvec1(2) posvec2(3) posvec1(4)]);

hold(ax{2}, 'on');
axis([0 inf 0.5 1.5]);

box(ax{2}, 'on');
% Set the remaining ax properties
set(ax{2}, 'XMinorTick', 'on', 'YMinorTick', 'on');

end
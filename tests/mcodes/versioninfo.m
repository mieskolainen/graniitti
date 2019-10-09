% Create version text

function h = versioninfo(jsonfile)

js = jsondecode(fileread(jsonfile));

xl = xlim; yl = ylim;
h{1} = text(xl(2)*1.25, yl(1)*1.025, sprintf('GRANIITTI %0.2f', js.version));
h{2} = text(xl(2)*1.25, yl(2)*0.060, 'github.com/mieskolainen');
for i = 1:2
set(h{i}, 'FontSize', 8, 'Color', ones(3,1) * 0.75, 'Rotation', 90);
end

end
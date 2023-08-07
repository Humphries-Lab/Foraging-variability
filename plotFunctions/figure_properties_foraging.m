% figure properties for foraging paper

% format = 'png'; % for panels tricky for EPS (e.g. Pcolor plots)
% color = 'rgb';
% dpi = 600;
% brewermap('demo') 

fontsize = 10;
fontname = 'Arial';
M = 3; % marker size for univariate scatter plots
sym = 'o';  % markers for scatters and strip plots

Units = 'centimeters';

% line widths
widths.plot = 1.5;
widths.error = 0.5;
widths.axis = 0.5;

% panel sizes
figsize.square = [20 20 7 7];

% % colours for matrix plots
% colourmaps.rate_map = brewermap(15,'Greys');
% colourmaps.inhibit_map = brewermap(15,'Blues');
% colourmaps.excite_map = brewermap(15,'OrRd');
% colourmaps.PC_map = brewermap(10,'PuOr');

% colours for lines
tmp = brewermap(9,'Set1');
color.rich = tmp(3,:);
color.poor = tmp(5,:);
color.patch = brewermap(3,'Blues');

color.text = [0.3 0.3 0.3];

% line types 
lines.model = '.--';
lines.exp = '.-';

% exportpath
export_path = '../outputs/';

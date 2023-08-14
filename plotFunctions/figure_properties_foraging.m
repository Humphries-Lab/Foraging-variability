% figure properties for foraging paper

% format = 'png'; % for panels tricky for EPS (e.g. Pcolor plots)
% color = 'rgb';
% dpi = 600;
% brewermap('demo') 

fontsize = 11;
fontname = 'Arial';
M_alpha = 0.6; % marker alpha for scatter plots

Units = 'centimeters';

% line widths
widths.plot = 2;
widths.error = 0.5;
widths.axis = 0.5;
widths.betaLines = 2;

% panel sizes
figsize.square = [20 20 7 7];
figsize.rectangle = [20 20 14 7];

% colours for lines
tmp_set1 = brewermap(9,'Set1');
color.rich = tmp_set1(2,:); % green
color.poor = tmp_set1(3,:); % purple
color.patch = brewermap(3,'Oranges');
color.general = tmp_set1(8,:); % pink

tmp_accent = brewermap(8,'Accent');
color.highlight = tmp_accent(6,:);
color.text = [0.3 0.3 0.3];

% line types 
lines.model = '.--';
lines.exp = '.-';

marker.rich = color.rich + (1-color.rich) * 0.2;
marker.poor = color.poor + (1-color.poor) * 0.2;

% exportpath
export_path = '../outputs/';

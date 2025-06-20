% figure properties for foraging paper

% format = 'png'; % for panels tricky for EPS (e.g. Pcolor plots)
% color = 'rgb';
% dpi = 600;
%brewermap('demo') 

fontsize = 11;
fontname = 'Arial';

Units = 'centimeters';

% line widths
widths.plot = 1.5;
widths.error = 0.5;
widths.axis = 0.5;
widths.betaLines = 2;

% panel sizes
figsize.small_panel = [20 20 4.5 4.5];
figsize.square = [20 20 6 6];
figsize.vertical = [20 20 3 6];
figsize.rectangle = [20 20 6 5];

% colours for lines
tmp_set1 = brewermap(8,'Set1');
color.poor =  tmp_set1(3,:); %green '#3AAA36'; 
color.rich =  tmp_set1(2,:); % blue '#2582C4'; 
tmp_oranges = brewermap(13,'Oranges');
color.patch = tmp_oranges([4, 7, 11],:);

color.general = '#F5A0B3'; % pink
color.highlight = '#9D004D';
color.scatter = '#857699';

color.text = [0.3 0.3 0.3];

% line types 
lines.model = '--';
lines.exp = '-';
lines.mvt = '-.';

% exportpath
export_path = '../';
overleaf_path = '/Users/exs165/Dropbox/Apps/Overleaf/240521 - foraging variability paper/figs/';

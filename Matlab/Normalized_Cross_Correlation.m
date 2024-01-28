clear all
close all
clc

%% read data
folder = 'Data';
rmpath(folder)
folder_data = 'Datanew';
addpath(folder_data)
data = readmatrix('PA0.csv');
data1 = readmatrix('PA6-18.csv');
data2 = readmatrix('PA6-18V.csv');
data3 = readmatrix('PA9-15.csv');
data4 = readmatrix('PA9-15V.csv');
data5 = readmatrix('PA9-27.csv');
data6 = readmatrix('PA9-27V.csv');

c{1}=data1;
c{2}=data2;
c{3}=data3;
c{4}=data4;
c{5}=data5;
c{6}=data6;

%% outline definition

x0 = data(:,1);
y0 = data(:,2);
k = boundary(x0, y0);   %obtains the boundary from the coordinates of the circle
writematrix ( [ x0(k), y0(k)] , 'circle_outline.csv'); %create the outline matrix

% trial = abs( data(:,5));   %questo è il dato che mi prende i valori a una certa f
% figure;
% scatter3( x0 , y0, trial );

% figure; 
% scatter ( x0 , y0 )
% hold on
% plot ( x0(k), y0(k), '-r', 'linewidth', 2);
% hold off



%% creo la mesh caso 0

gridPointsX = 80;
gridPointsY = 80;

x_max0 = max(x0); x_min0 = min(x0); y_max0 = max(y0); y_min0 = min(y0);
xx0 = linspace ( x_min0 , x_max0 , gridPointsX );
yy0 = linspace ( y_min0, y_max0, gridPointsY );
[X0, Y0] = meshgrid (xx0, yy0);

% figure;
% scatter ( x0 , y0 ) ; hold on; scatter ( X0(:) , Y0(:) , '*r');

%% taglio ciò che è fuori dal cerchio
outline0 = readmatrix ( 'circle_outline.csv');

% define the inner circle
in0 = inpolygon( X0(:), Y0(:), outline0(:,1), outline0(:,2));
IN0 = reshape( in0 , [gridPointsY, gridPointsX]);   % create binary mask inside the boundary

%% NCC

NCCmatr = zeros(499, 6);

for i=1:499   %take the freq of normal cajos and us it to confront with the others
    trial = abs(data(:,i));
    trial (isnan(trial)) = 0;
    Vq0 = griddata ( x0 , y0 , trial , X0 , Y0 , 'cubic');
    press = Vq0.*IN0;
    press(isnan(press)) = 0; 
    press( press == 0 ) = [];
    I = press(:);
    for j = 1:6
        Case = c{j};
        x = Case(:,1); 
        y = Case(:,2);        
        t = abs(Case(:,i));
        t (isnan(t)) = 0;
        Vq = griddata ( x , y , t , X0 , Y0 , 'cubic');
        pression = Vq.*IN0;
        pression(isnan(pression)) = 0; 
        pression( pression == 0 ) = [];
        G = pression(:);
        NCCmatr(i,j) = (abs(I)'*abs(G))/(norm(I)*norm(G));
    end
end

%% plots


figure(1);
imagesc(NCCmatr')
h = colorbar;
ylabel(h, '|c|')
set(h.Label,'FontSize',15)
%title('NCC matrix')
yticklabels({'6-18', '6-18 V', '9-15', '9-15 V', '9-27', '9-27 V'})
xticks([1:50:5000])
% xticklabels({'', '', '', '', '', '', '', '', '100', '', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '500','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '1000','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '1500','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '','','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '2000','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '2500','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '3000','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '3500','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '4000','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '4500','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '','', '', '', '', '5000'})
xticklabels({'20', '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500', '5000'})
ax = gca;
ax.FontSize = 14; 
ylabel('Metamaterial case', 'FontSize', 14)
xlabel('Frequency [Hz]', 'FontSize', 14)
xtickangle( ax , 90 )
pbaspect([1.25 1 1])
% exportgraphics(gcf, 'test.pdf'); % use this to save as a cropped pdf

clear all;
close all;

%%
%rifaccio il cut plane spostato dalla superficie

folder_data = 'Data';
addpath(folder_data)
data = readmatrix('PA0.csv');
x = data(:,1);
y = data(:,2);
trial = abs( data(:,3));   %questo è il dato che mi prende i valori a una certa f
% P = inpaint_nans(trial);

figure;
scatter3( x , y, trial );

idx = find ( isnan(trial ) == 1 ) ;
 if isempty ( idx ) ~= 1
    figure;
    scatter ( x , y)
    hold on;
    scatter ( x (idx) , y( idx ))
 end
 
%% creo la mesh
gridPointsX = 80;
gridPointsY = 80;

x_max = max(x); x_min = min(x); y_max = max(y); y_min = min(y);
xx = linspace ( x_min , x_max , gridPointsX );
yy = linspace ( y_min, y_max, gridPointsY );
[X, Y] = meshgrid (xx, yy);

figure;
scatter ( x , y ) ; hold on; scatter ( X(:) , Y(:) , '*r');


%% interpolation
trial (isnan(trial)) = 0;
Vq = griddata ( x , y , trial , X , Y , 'cubic');

figure; 
mesh ( X , Y , Vq ) ; hold on; scatter3( x , y , trial , '*' ); 


%% taglio fuori dal cerchio
outline = readmatrix ( 'circle_outline.csv');

% define the inner polygon
in = inpolygon( X(:), Y(:), outline(:,1), outline(:,2));
IN = reshape( in , [gridPointsY, gridPointsX]);   %creo la maschera binaria solo dentro al cerchio, vale per tutti

press = Vq.*IN; %ottengo i valori solo dentro al cerchio, fuori sono tutti 0

press ( press == 0 ) = [];  %elimino tutti gli zeri fuori dal cerchio
press = press(:);   %questo è il dato che vado ad usare per la NCC


%% NCC

NCCmatr = zeros(49, 6);

for i=1:49   %prendo la frequenza dal caso 0 e lo uso per confrontare gli altri
    trial = abs(data(:,i+2));
    trial (isnan(trial)) = 0;
    Vq = griddata ( x , y , trial , X , Y , 'cubic');
    press = Vq.*IN;
    press( press == 0 ) = [];
    I = press(:);
    for j = 1:6
        Y = c{j};
        Y = Y(:,i);
        NCCmatr(i,j) = (abs(I)'*abs(Y))/(norm(I)*norm(Y));
    end
end

%% plots

figure(1);
imagesc(NCCmatr)
colorbar
title('Correlation matrix')
xticklabels({'6-18', '6-18 V', '9-15', '9-15 V', '9-27', '9-27 V'})
yticks([1:1:49])
yticklabels({'20', '22.4', '25', '28', '31.5', '35.5', '40', '45', '50', '56', '63', '71', '80', '90', '100', '112', '125', '140', '160', '180', '200', '224', '250', '280', '315', '355', '400', '450', '500', '560', '630', '710', '800', '900', '1000', '1120', '1250', '1400', '1600', '1800', '2000', '2240', '2500', '2800', '3150', '3550', '4000', '4500', '5000'})
xlabel('Metamaterial case')
ylabel('Frequency [Hz]')



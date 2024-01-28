data = readmatrix ( 'PA0.csv');
x = data(:,1);
y = data(:,2);
k = boundary(x, y);   %ottiene il confine del cerchio dalle coordinate
writematrix ( [ x(k), y(k)] , 'circle_outline.csv'); %crea la matrice dell'outline

figure; 
scatter ( x , y )
hold on
plot ( x(k), y(k), '-r', 'linewidth', 2);
hold off


%%
% %% 6-18
% data1 = readmatrix ( 'PA6-18.csv');
% x1 = data1(:,1);
% y1 = data1(:,2);
% k1 = boundary(x1, y1);  
% k1 = interp1((1:numel(k1)), k1, linspace(1, numel(k1), numel(k)), 'nearest')';
% writematrix ( [ x1(k1), y1(k1)] , 'circle_outline1.csv');
% %% 6-18V
% data2 = readmatrix ( 'PA6-18V.csv');
% x2 = data2(:,1);
% y2 = data2(:,2);
% k2= boundary(x2, y2); 
% k2 = interp1((1:numel(k2)), k2, linspace(1, numel(k2), numel(k)), 'nearest')';
% writematrix ( [ x2(k2), y2(k2)] , 'circle_outline2.csv');
% %% 9-15
% data3 = readmatrix ( 'PA9-15.csv');
% x3 = data3(:,1);
% y3 = data3(:,2);
% k3 = boundary(x3, y3);  
% k3 = interp1((1:numel(k3)), k3, linspace(1, numel(k3), numel(k)), 'nearest')';
% writematrix ( [ x3(k3), y3(k3)] , 'circle_outline3.csv');
% %% 9-15V
% data4 = readmatrix ( 'PA9-15V.csv');
% x4 = data3(:,1);
% y4= data3(:,2);
% k4 = boundary(x4, y4); 
% k4 = interp1((1:numel(k4)), k4, linspace(1, numel(k4), numel(k)), 'nearest')';
% writematrix ( [ x4(k4), y4(k4)] , 'circle_outline4.csv');
% %% 9-27
% data5 = readmatrix ( 'PA9-27.csv');
% x5 = data3(:,1);
% y5= data3(:,2);
% k5 = boundary(x5, y5);   
% k5 = interp1((1:numel(k5)), k5, linspace(1, numel(k5), numel(k)), 'nearest')';
% writematrix ( [ x5(k5), y5(k5)] , 'circle_outline5.csv');
% %% 9-27V
% data6 = readmatrix ( 'PA9-27V.csv');
% x6= data3(:,1);
% y6 = data3(:,2);
% k6 = boundary(x6, y6);   
% k6 = interp1((1:numel(k6)), k6, linspace(1, numel(k6), numel(k)), 'nearest')';
% writematrix ( [ x6(k6), y6(k6)] , 'circle_outline6.csv');
% 
% %% interpolazione 
% %k è 272 nel caso 0, 271 nel caso 1 e 270 in tutti gli altri
% %non posso appendere 0 perchè non va bene a write matrix --> interpolo (con nearest come metodo)

% d{1}=readmatrix ('circle_outline1.csv');
% d{2}=readmatrix ('circle_outline2.csv');
% d{3}=readmatrix ('circle_outline3.csv');
% d{4}=readmatrix ('circle_outline4.csv');
% d{5}=readmatrix ('circle_outline5.csv');
% d{6}=readmatrix ('circle_outline6.csv');




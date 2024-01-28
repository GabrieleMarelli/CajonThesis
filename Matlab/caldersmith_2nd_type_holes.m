close all
clear all
clc

%use matlab from R2020a version
%CHANGE THE MATRIX F USED FOR HORIZONTAL OR VERTICAL
%% data

rho = 620; %density plywood
rhoA = 1.225; %density air
a = 0.425;
b = 0.299;
h = 0.003;
t = 0.0025;    %thickness hole
N = 50;     %number of holes

R = [0.0045, 0.0075     %semi axis y and x ellipse  3-5 (9mm, 15mm)
     0.003, 0.009       %3-9 (6mm, 18mm)
     0.0045, 0.0135];   %3-9 (9mm, 27mm)
 
dim = size(R,1);

%% Calculations

volM = zeros(dim,1);
volDiff = zeros(dim,1);
massaDiff = zeros(dim,1);
massaM = zeros(dim,1);
massaTot = zeros(dim,1);
rhoh = zeros(dim,1);

volTot = a*b*h;                         %total volume of the plate
for k = 1:dim
    volM(k) = pi*R(k,1)*R(k,2)*t*N;     %total volume of the holes
    volDiff(k) = volTot - volM(k);      %volume only plywood
    massaDiff(k) = volDiff(k)*rho;      %mass only plywood
    massaM(k) = volM(k)*rhoA;           %mass only air
    massaTot(k) = massaDiff(k) + massaM(k);
    rhoh(k) = massaTot(k)/volTot;
end

%% natural freq of the tapa, from comsol eigenfreq
%    f02,    f20,     f11
% F = [72.304, 106.15, 30.542        %normal tapa
%      70.13, 104.32, 30.225         %tapa 9mm-15mm 2ND TYPE holes HORIZONTAL
%      69.29, 105.09, 30.219         %tapa 6mm-18mm 2ND TYPE HORIZONTAL
%      67.087, 104.29, 30.078];      %tapa 9mm-27mm 2ND TYPE HORIZONTAL
 
F = [72.304, 106.15, 30.542       %normal tapa
     71.423, 102.38, 30.243       %tapa 9mm-15mm 2ND TYPE holes VERTICAL
     71.877, 101.19, 30.26        %tapa 6mm-18mm 2ND TYPE VERTICAL
     71.61,  95.002, 29.93];      %tapa 9mm-27mm 2ND TYPE VERTICAL

%% Caldersmith's formulae  (constants for elastic deformation)
Dh = zeros(3,dim);
Eh = zeros(2,dim);
Gh = zeros(1,dim);

%true values
El = 11140*10^6;
Er = 5880*10^6;
Glr = 910*10^6;

%values for the normal tapa
D = [0.08006, 0.08006, 0.274]'.*[a^4, b^4, a^2*b^2]'.*(F(1,:)'.^2)*rho./(h^2); 
E = [12*D(1,1); 12*D(2,1)];
G = 3*D(3,1);

%error calculation
El_err = (E(1,1)-El)/El*100;
Er_err = (E(2,1)-Er)/Er*100;
fprintf('Error bars: El = %.1f%%,  Er = %.1f%%\n', El_err, Er_err);
Glr_err = abs((G-Glr)/Glr)*100;
fprintf('Error bars: Glr = %.1f%%\n\n', Glr_err);

%values for the tapa with holes
for j = 1:dim
    Dh(:,j) = [0.08006, 0.08006, 0.274]'.*[a^4, b^4, a^2*b^2]'.*(F(j+1,:)'.^2)*rhoh(j)./(h^2); 
    Eh(1,j) = 12*Dh(1,j);
    Eh(2,j) = 12*Dh(2,j);
    Gh(1,j) = 3*Dh(3,j);
end


diff = zeros(3,1);
diffE = zeros(2,1);
index = [1,3,4];
for i = 1:dim
    fprintf('case %d: holes of %.2fmm-%.2fmm\n',i,2*R(i,1)*10^3,2*R(i,2)*10^3)
    fprintf('New density %d = %.1f\nDensity variation = %.1f%%\n',i, rhoh(i), (1- rhoh(i)/rho)*100);
    for m = 1:3
        diff(m,i) = (1- Dh(m,i)/D(m,1))*100;
        fprintf('D%d = %d, D%dh = %d\n',index(m),D(m,1),index(m),Dh(m,i))
        fprintf('Percentage variation D%d = %.1f%%\n',index(m),diff(m,i))
    end
    fprintf('E1 = %d, E1h = %d\nE2 = %d, E2h = %d\n',E(1,1),Eh(1,i),E(2,1),Eh(2,i))
    for n = 1:2
        diffE(n,i) = (1- Eh(n,i)/E(n,1))*100;
        fprintf('Percentage variation E%d = %.1f%%\n',n,diffE(n,i))
    end
    diffG = (1- Gh(1,i)/G)*100;
    
    fprintf('G = %d, Gh = %d\n',G,Gh(1,i))
    fprintf('Percentage variation G = %.1f%%\n\n',diffG)
    
end

%% Plots

figure(1)
errorbar(rhoh(:,1), Eh(1,:),El_err*ones(size(rhoh))*El/100,'o', 'MarkerSize',8, 'MarkerFaceColor', 'blue');
% scatter(rhoh(:,1), Eh(1,:),'DisplayName', 'E1 holes 2nd type');
% hold on 
% errorbar(rhoh(:,1), Eh(2,:),Er_err*ones(size(rhoh))*Er/100,'o','DisplayName', 'E2 holes 2nd type');
% scatter(rhoh(:,1), Eh(2,:),'DisplayName', 'E2 holes 2nd type');
yline(El, 'b--', 'LineWidth', 1);
% yline(Er, 'r--', 'LineWidth', 0.5,'DisplayName', 'E2 normal tapa');
lgd = legend('E1 holes', 'E1 normal tapa');
lgd.FontSize = 12;
lgd.Title.String = '2nd type holes'; 
set(lgd, 'Box', 'off');
ax = gca;
ax.FontSize = 12; 
xlabel('Density [kg/m^3]', 'FontSize', 13)
ylabel('Young''s modulus E [Pa]', 'FontSize', 13)
% hold off


figure(2)
e = errorbar(rhoh(:,1), Eh(2,:),Er_err*ones(size(rhoh))*Er/100,'o', 'MarkerSize',8,'MarkerFaceColor', 'red');
e.Color = 'red';
yline(Er, 'r--', 'LineWidth', 1);
lgd = legend('E2 holes', 'E2 normal tapa');
lgd.FontSize = 12;
lgd.Title.String = '2nd type holes'; 
set(lgd, 'Box', 'off');
ax = gca;
ax.FontSize = 12; 
xlabel('Density [kg/m^3]', 'FontSize', 13)
ylabel('Young''s modulus E [Pa]', 'FontSize', 13)


figure(3)
e1 = errorbar(rhoh(:,1), Gh(1,:),Glr_err*ones(size(rhoh))*Glr/100,'o', 'MarkerSize',8,'MarkerFaceColor', 'magenta');
e1.Color = 'magenta';
% scatter(rhoh(:,1), Gh(1,:),'DisplayName', 'G holes 2nd type');
yline(Glr, 'm--', 'LineWidth', 1);
lgd = legend('G holes', 'G normal tapa');
lgd.FontSize = 12;
lgd.Title.String = '2nd type holes'; 
set(lgd, 'Box', 'off');
ax = gca;
ax.FontSize = 12; 
xlabel('Density [kg/m^3]', 'FontSize', 13)
ylabel('Shear modulus G [Pa]', 'FontSize', 13)
hold off

% figure(4)
% scatter(rhoh(:,1), Dh(1,:),'DisplayName', 'D1 holes 2nd type')
% hold on
% scatter(rhoh(:,1), Dh(2,:),'DisplayName', 'D3 holes 2nd type')
% hold on
% scatter(rhoh(:,1), Dh(3,:),'g','DisplayName', 'D4 holes 2nd type')
% yline(D(1,1), 'b--', 'LineWidth', 0.5,'DisplayName', 'D1 normal tapa');
% yline(D(2,1), 'r--', 'LineWidth', 0.5,'DisplayName', 'D3 normal tapa');
% yline(D(3,1), 'g--', 'LineWidth', 0.5,'DisplayName', 'D4 normal tapa');
% legend
% xlabel('density [kg/m^3]')
% ylabel('Elastic deformation D [Pa]')
% hold off





% figure(1)
% scatter(rhoh(:,1), Eh(1,:),'DisplayName', 'E1 2nd type holes');
% hold on 
% scatter(rhoh(:,1), Eh(2,:),'DisplayName', 'E2 2nd type holes');
% yline(E(1,1), 'b--', 'LineWidth', 0.5,'DisplayName', 'E1 normal tapa');   %elastic modulus longitudinal of normal tapa
% yline(E(2,1), 'r--', 'LineWidth', 0.5, 'DisplayName', 'E2 normal tapa');   %elastic modulus radial of normal tapa
% legend
% xlabel('density [kg/m^3]')
% ylabel('Young''s modulus [Pa]')
% hold off
% 
% figure(2)
% scatter(rhoh(:,1), Gh(1,:),'DisplayName', 'G modulus tapa with 2nd type holes')
% yline(G(1,1), 'b--', 'LineWidth', 0.5,'DisplayName', 'G modulus normal tapa');
% legend
% xlabel('density [kg/m^3]')
% ylabel('Shear modulus [Pa]')
% 
% figure(3)
% scatter(rhoh(:,1), Dh(1,:),'DisplayName', 'D1 2nd type holes')
% hold on
% scatter(rhoh(:,1), Dh(2,:),'DisplayName', 'D3 2nd type holes')
% hold on
% scatter(rhoh(:,1), Dh(3,:),'g','DisplayName', 'D4 2nd type holes')
% yline(D(1,1), 'b--', 'LineWidth', 0.5,'DisplayName', 'D1 normal tapa');
% yline(D(2,1), 'r--', 'LineWidth', 0.5,'DisplayName', 'D3 normal tapa');
% yline(D(3,1), 'g--', 'LineWidth', 0.5,'DisplayName', 'D4 normal tapa');
% legend
% xlabel('density [kg/m^3]')
% ylabel('Elastic deformation D [Pa]')
% hold off



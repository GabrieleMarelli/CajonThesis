close all
clear all
clc

%use matlab from R2020a version
%% data

rho = 620; %density plywood
a = 0.425;
b = 0.299;
h = 0.003;


%% natural freq of the tapa, from comsol eigenfreq
%    f02,    f20,     f11
F = [72.304, 106.15, 30.542];        %normal tapa
   

%% Caldersmith's formulae  (constants for elastic deformation)
%true values
El = 11140*10^6;
Er = 5880*10^6;
Glr = 910*10^6;

%values for the normal tapa
D = [0.08006, 0.08006, 0.274]'.*[a^4, b^4, a^2*b^2]'.*(F(1,:)'.^2)*rho./(h^2); 
E = [12*D(1,1); 12*D(2,1)];
G = 3*D(3,1);

%error calculation
El_err = abs((El - E(1,1)))/El*100;
Er_err = abs((Er - E(2,1)))/Er*100;
fprintf('Error bars: El = %.1f%%,  Er = %.1f%%\n', El_err, Er_err);
Glr_err = abs((Glr - G)/Glr)*100;
fprintf('Error bars: Glr = %.1f%%\n\n', Glr_err);


%% Plots

figure(1)
%scatter(rho, E(1,1),'DisplayName', 'E1  Caldersmith');
errorbar(rho, E(1,1),El_err*El/100,'o','DisplayName', 'E1 Caldersmith');
hold on
yline(El, 'b--', 'LineWidth', 0.5,'DisplayName', 'E1 REAL normal tapa');
legend
xlabel('density [kg/m^3]')
ylabel('Young''s modulus E [Pa]')
hold off

figure(2)
% scatter(rho, E(2,1),'DisplayName', 'E2 Caldersmith');
errorbar(rho, E(2,1),Er_err*Er/100,'o','DisplayName', 'E2 Caldersmith');
hold on
yline(Er, 'b--', 'LineWidth', 0.5,'DisplayName', 'E2 REAL normal tapa');
legend
xlabel('density [kg/m^3]')
ylabel('Young''s modulus E [Pa]')
hold off

figure(3)
%scatter(rhoh(:,1), Gh(1,:),'DisplayName', 'G modulus tapa with holes')
errorbar(rho, G(1,1),Glr_err/100*Glr,'o','DisplayName', 'G tapa Caldersmith');
yline(Glr, 'b--', 'LineWidth', 0.5,'DisplayName', 'G modulus REAL normal tapa');
legend
xlabel('density [kg/m^3]')
ylabel('Shear modulus [Pa]')
ylim([8.0*10^8 9.2*10^8])
hold off


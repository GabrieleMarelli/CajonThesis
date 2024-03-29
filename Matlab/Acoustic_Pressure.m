close all
clear all
clc
%%

% data0 = [20.000	2.0897E-5-5.7107E-7i
% 22.400	2.1106E-5-1.7566E-6i
% 25.000	1.8398E-5-2.3916E-6i
% 28.000	1.0070E-5-1.6398E-6i
% 31.500	-9.6592E-6+2.5098E-6i
% 35.500	-5.1913E-5+1.4390E-5i
% 40.000	-1.3910E-4+4.3629E-5i
% 45.000	-3.1926E-4+1.1307E-4i
% 50.000	-6.5987E-4+2.6070E-4i
% 56.000	-0.0015310+6.8308E-4i
% 63.000	-0.0046805+0.0023917i
% 71.000	-0.36964+0.70510i
% 80.000	0.010983-0.0072798i
% 90.000	0.0081826-0.0063963i
% 100.00	0.0081911-0.0075257i
% 112.00	0.0097187-0.011055i
% 125.00	0.015025-0.025364i
% 140.00	-0.023840+0.020437i
% 160.00	-1.2759E-4+1.7430E-4i
% 180.00	8.5335E-4-0.0019001i
% 200.00	0.0017092-0.0065888i
% 224.00	-6.2220E-4+0.038093i
% 250.00	-0.0020247+0.026048i
% 280.00	-0.0036091-0.018042i
% 315.00	-0.032434-0.073472i
% 355.00	0.0069554+0.0097222i
% 400.00	0.010115+0.017871i
% 450.00	0.66676+1.3689i
% 500.00	0.0099451+0.0088767i
% 560.00	-0.024117-0.010998i
% 630.00	0.0012491-0.0032728i
% 710.00	-0.0055541-0.0031688i
% 800.00	0.011510-0.019295i
% 900.00	-0.0011560-0.0076158i
% 1000.0	-0.085736-0.049059i
% 1120.0	0.12215+0.071009i
% 1250.0	-0.0092968+0.030321i
% 1400.0	-0.087820+0.021121i
% 1600.0	0.0099074+0.078740i
% 1800.0	0.0028226-0.0067367i
% 2000.0	0.032942-0.020155i
% 2240.0	0.0012369+0.010951i
% 2500.0	0.031642-0.031251i
% 2800.0	0.035755+0.0020176i
% 3150.0	0.047799-2.4353E-4i
% 3550.0	0.24049-0.17726i
% 4000.0	0.25904+0.051084i
% 4500.0	0.0069647-0.0040399i
% 5000.0	-0.011216+0.039324i
% ];
% 
% %9-15
% data1 = [20.000	2.0725E-5-5.6061E-7i
% 22.400	2.0804E-5-1.7227E-6i
% 25.000	1.7892E-5-2.3121E-6i
% 28.000	9.2252E-6-1.4715E-6i
% 31.500	-1.1082E-5+2.8498E-6i
% 35.500	-5.4320E-5+1.5056E-5i
% 40.000	-1.4316E-4+4.4908E-5i
% 45.000	-3.2602E-4+1.1547E-4i
% 50.000	-6.7017E-4+2.6478E-4i
% 56.000	-0.0015439+6.8883E-4i
% 63.000	-0.0046418+0.0023716i
% 71.000	-0.16251+0.12769i
% 80.000	0.011591-0.0076726i
% 90.000	0.0084813-0.0066231i
% 100.00	0.0084284-0.0077324i
% 112.00	0.0099122-0.011233i
% 125.00	0.014988-0.024635i
% 140.00	-0.030720+0.024129i
% 160.00	-5.0596E-4+8.0594E-4i
% 180.00	5.2555E-4-0.0010529i
% 200.00	0.0016381-0.0062797i
% 224.00	-6.4056E-4+0.038515i
% 250.00	-0.0020349+0.026057i
% 280.00	-0.0035110-0.017664i
% 315.00	-0.028065-0.064181i
% 355.00	0.0084456+0.011880i
% 400.00	0.012065+0.020332i
% 450.00	-0.97154-0.14825i
% 500.00	0.033588+0.012568i
% 560.00	-0.013133-0.010780i
% 630.00	0.0027234-0.0037903i
% 710.00	-0.0029672-0.0049208i
% 800.00	0.015099-0.023243i
% 900.00	-0.0015365-0.0083715i
% 1000.0	-0.058390-0.035555i
% 1120.0	0.13857+0.071166i
% 1250.0	-0.0036810+0.014269i
% 1400.0	-0.089777+0.019693i
% 1600.0	0.016601+0.10917i
% 1800.0	-0.031263-0.0014874i
% 2000.0	0.022014-0.011235i
% 2240.0	-0.012598+0.014926i
% 2500.0	0.038185-0.084826i
% 2800.0	0.035202-0.0044224i
% 3150.0	0.046816-0.0029211i
% 3550.0	0.24212-0.20745i
% 4000.0	-0.0018179-0.44653i
% 4500.0	-0.021202-0.031242i
% 5000.0	-0.0050861+0.028730i
% ];
% 
% %9-15V
% data2 = [20.000	2.0486E-5-5.5209E-7i
% 22.400	2.0496E-5-1.6935E-6i
% 25.000	1.7497E-5-2.2545E-6i
% 28.000	8.7105E-6-1.3724E-6i
% 31.500	-1.1771E-5+3.0124E-6i
% 35.500	-5.5269E-5+1.5318E-5i
% 40.000	-1.4450E-4+4.5334E-5i
% 45.000	-3.2799E-4+1.1617E-4i
% 50.000	-6.7304E-4+2.6593E-4i
% 56.000	-0.0015483+6.9079E-4i
% 63.000	-0.0046461+0.0023737i
% 71.000	-0.15435+0.11909i
% 80.000	0.011654-0.0077147i
% 90.000	0.0085142-0.0066490i
% 100.00	0.0084567-0.0077584i
% 112.00	0.0099420-0.011266i
% 125.00	0.015020-0.024661i
% 140.00	-0.031374+0.024474i
% 160.00	-5.8681E-4+9.4148E-4i
% 180.00	5.3318E-4-0.0010695i
% 200.00	0.0016395-0.0062836i
% 224.00	-6.4379E-4+0.037963i
% 250.00	-0.0020168+0.025841i
% 280.00	-0.0034738-0.017485i
% 315.00	-0.027491-0.062969i
% 355.00	0.0085915+0.012102i
% 400.00	0.012836+0.021732i
% 450.00	-1.0575-0.12890i
% 500.00	0.030757+0.012023i
% 560.00	-0.012888-0.010862i
% 630.00	0.0026951-0.0038087i
% 710.00	-0.0027887-0.0049225i
% 800.00	0.014880-0.024211i
% 900.00	-9.6257E-4-0.0067240i
% 1000.0	-0.056622-0.030986i
% 1120.0	0.12973+0.068546i
% 1250.0	-0.0014172+0.0093014i
% 1400.0	-0.089927+0.019029i
% 1600.0	0.015227+0.10967i
% 1800.0	-0.040917-0.0072774i
% 2000.0	0.021944-0.014695i
% 2240.0	0.0018804-0.0039777i
% 2500.0	0.043363-0.071284i
% 2800.0	0.038007-0.0041061i
% 3150.0	0.045760-0.0015856i
% 3550.0	0.25897-0.13939i
% 4000.0	0.035996-0.28975i
% 4500.0	-0.028945-0.033028i
% 5000.0	-0.0056042+0.026176i];
% 
% %6-18
% data3 = [20.000	2.0592E-5-5.5868E-7i
% 22.400	2.0673E-5-1.7136E-6i
% 25.000	1.7784E-5-2.3000E-6i
% 28.000	9.1855E-6-1.4673E-6i
% 31.500	-1.0965E-5+2.8194E-6i
% 35.500	-5.3863E-5+1.4929E-5i
% 40.000	-1.4198E-4+4.4539E-5i
% 45.000	-3.2329E-4+1.1450E-4i
% 50.000	-6.6431E-4+2.6247E-4i
% 56.000	-0.0015290+6.8218E-4i
% 63.000	-0.0045847+0.0023422i
% 71.000	-0.14381+0.10887i
% 80.000	0.011577-0.0076641i
% 90.000	0.0084431-0.0065941i
% 100.00	0.0083795-0.0076889i
% 112.00	0.0098435-0.011160i
% 125.00	0.014868-0.024490i
% 140.00	-0.029716+0.023490i
% 160.00	-3.5340E-4+5.4568E-4i
% 180.00	5.3644E-4-0.0010816i
% 200.00	0.0016479-0.0063197i
% 224.00	-6.4114E-4+0.038338i
% 250.00	-0.0020457+0.026158i
% 280.00	-0.0035102-0.017659i
% 315.00	-0.028920-0.066046i
% 355.00	0.0080496+0.011305i
% 400.00	0.011753+0.019904i
% 450.00	-0.99143-0.14588i
% 500.00	0.032057+0.012287i
% 560.00	-0.013249-0.010535i
% 630.00	0.0023362-0.0037539i
% 710.00	-0.0034159-0.0046137i
% 800.00	0.014499-0.023415i
% 900.00	-9.5178E-4-0.0071732i
% 1000.0	-0.062434-0.037522i
% 1120.0	0.13224+0.069689i
% 1250.0	-0.0042336+0.015242i
% 1400.0	-0.083313+0.015123i
% 1600.0	0.014902+0.10207i
% 1800.0	-0.026123-9.8119E-4i
% 2000.0	0.030544-0.013834i
% 2240.0	-0.0056317+0.0091398i
% 2500.0	0.036786-0.047450i
% 2800.0	0.030648-0.0025524i
% 3150.0	0.044733-0.0015960i
% 3550.0	0.25975-0.19657i
% 4000.0	0.57979-0.20629i
% 4500.0	-0.021072-0.028670i
% 5000.0	-0.0052586+0.035503i
% ];
% 
% %6-18V
% data4 = [20.000	2.0808E-5-5.6491E-7i
% 22.400	2.0935E-5-1.7367E-6i
% 25.000	1.8095E-5-2.3435E-6i
% 28.000	9.5486E-6-1.5354E-6i
% 31.500	-1.0555E-5+2.7243E-6i
% 35.500	-5.3442E-5+1.4813E-5i
% 40.000	-1.4167E-4+4.4439E-5i
% 45.000	-3.2344E-4+1.1455E-4i
% 50.000	-6.6587E-4+2.6308E-4i
% 56.000	-0.0015363+6.8543E-4i
% 63.000	-0.0046314+0.0023663i
% 71.000	-0.17904+0.14556i
% 80.000	0.011469-0.0075932i
% 90.000	0.0084163-0.0065735i
% 100.00	0.0083758-0.0076861i
% 112.00	0.0098705-0.011194i
% 125.00	0.015004-0.024788i
% 140.00	-0.029163+0.023388i
% 160.00	-4.8210E-4+7.7057E-4i
% 180.00	5.9301E-4-0.0012251i
% 200.00	0.0016441-0.0063084i
% 224.00	-6.4287E-4+0.038655i
% 250.00	-0.0020390+0.026145i
% 280.00	-0.0034902-0.017531i
% 315.00	-0.027659-0.063274i
% 355.00	0.0083307+0.011712i
% 400.00	0.012001+0.020490i
% 450.00	-1.8752+0.56536i
% 500.00	0.026760+0.011241i
% 560.00	-0.015218-0.010954i
% 630.00	0.0023818-0.0036895i
% 710.00	-0.0033329-0.0046442i
% 800.00	0.014679-0.023208i
% 900.00	-0.0010162-0.0070585i
% 1000.0	-0.057186-0.033710i
% 1120.0	0.12786+0.067962i
% 1250.0	-0.0028889+0.012298i
% 1400.0	-0.089100+0.020250i
% 1600.0	0.012476+0.10207i
% 1800.0	-0.029189-0.011200i
% 2000.0	0.035032-0.015155i
% 2240.0	-0.012290+0.0024282i
% 2500.0	0.037955-0.041617i
% 2800.0	0.036355+7.9260E-4i
% 3150.0	0.046404-0.0012328i
% 3550.0	0.27156-0.17341i
% 4000.0	0.32323-0.44586i
% 4500.0	-0.011072-0.050850i
% 5000.0	-0.016621+0.023836i
% ];
% 
% %9-27
% data5 = [20.000	2.0634E-5-5.5631E-7i
% 22.400	2.0632E-5-1.7047E-6i
% 25.000	1.7593E-5-2.2662E-6i
% 28.000	8.7170E-6-1.3711E-6i
% 31.500	-1.1945E-5+3.0556E-6i
% 35.500	-5.5773E-5+1.5458E-5i
% 40.000	-1.4556E-4+4.5665E-5i
% 45.000	-3.2982E-4+1.1682E-4i
% 50.000	-6.7533E-4+2.6683E-4i
% 56.000	-0.0015472+6.9029E-4i
% 63.000	-0.0045888+0.0023443i
% 71.000	-0.10311+0.072360i
% 80.000	0.012088-0.0079924i
% 90.000	0.0086999-0.0067882i
% 100.00	0.0085862-0.0078678i
% 112.00	0.010013-0.011314i
% 125.00	0.014843-0.023929i
% 140.00	-0.038315+0.027243i
% 160.00	-8.1713E-4+0.0013186i
% 180.00	3.5641E-4-6.1573E-4i
% 200.00	0.0016278-0.0062257i
% 224.00	-6.4315E-4+0.038201i
% 250.00	-0.0020234+0.025657i
% 280.00	-0.0035078-0.017770i
% 315.00	-0.026722-0.061382i
% 355.00	0.0097138+0.013739i
% 400.00	0.012454+0.020475i
% 450.00	-0.40818-0.15445i
% 500.00	0.051677+0.016922i
% 560.00	-0.0073073-0.011007i
% 630.00	0.0037203-0.0041840i
% 710.00	-0.0014113-0.0061575i
% 800.00	0.017280-0.025770i
% 900.00	-0.012882-0.055515i
% 1000.0	-0.052656-0.033333i
% 1120.0	0.15703+0.072891i
% 1250.0	-2.6709E-4+0.0055189i
% 1400.0	-0.085598+0.014126i
% 1600.0	0.033391+0.14551i
% 1800.0	-0.063116+0.0092671i
% 2000.0	0.029051-0.0090505i
% 2240.0	-0.012783+0.031829i
% 2500.0	-0.031258+0.014682i
% 2800.0	-0.0046731-0.034558i
% 3150.0	0.049359-0.0056420i
% 3550.0	0.13447-0.22594i
% 4000.0	-0.061343-0.25673i
% 4500.0	-0.045899-0.13427i
% 5000.0	-0.037311+0.019308i];
% 
% %9-27V
% data6 = [20.000	2.0450E-5-5.4580E-7i
% 22.400	2.0363E-5-1.6750E-6i
% 25.000	1.7199E-5-2.2045E-6i
% 28.000	8.1181E-6-1.2520E-6i
% 31.500	-1.2897E-5+3.2826E-6i
% 35.500	-5.7357E-5+1.5897E-5i
% 40.000	-1.4833E-4+4.6539E-5i
% 45.000	-3.3499E-4+1.1866E-4i
% 50.000	-6.8521E-4+2.7075E-4i
% 56.000	-0.0015706+7.0078E-4i
% 63.000	-0.0046788+0.0023905i
% 71.000	-0.12315+0.089541i
% 80.000	0.012039-0.0079625i
% 90.000	0.0087215-0.0068061i
% 100.00	0.0086290-0.0079082i
% 112.00	0.010088-0.011401i
% 125.00	0.015009-0.024223i
% 140.00	-0.038591+0.027570i
% 160.00	-9.1831E-4+0.0014882i
% 180.00	1.9663E-4-2.0467E-4i
% 200.00	0.0015674-0.0059775i
% 224.00	-6.6961E-4+0.039385i
% 250.00	-0.0020125+0.025603i
% 280.00	-0.0034410-0.017414i
% 315.00	-0.024496-0.056488i
% 355.00	0.010135+0.014352i
% 400.00	0.014124+0.023372i
% 450.00	-1.0513-0.11975i
% 500.00	0.041924+0.014505i
% 560.00	-0.0089299-0.011583i
% 630.00	0.0038675-0.0041382i
% 710.00	-1.7639E-4-0.0066292i
% 800.00	0.018930-0.027534i
% 900.00	-0.0028869-0.011529i
% 1000.0	-0.039402-0.023831i
% 1120.0	0.15718+0.073325i
% 1250.0	0.0066216-0.010226i
% 1400.0	-0.098990+0.021531i
% 1600.0	0.028482+0.14836i
% 1800.0	-0.048502-0.079108i
% 2000.0	0.023838-0.0064689i
% 2240.0	9.0128E-4+0.027867i
% 2500.0	-0.11306+0.030938i
% 2800.0	0.097280-0.027743i
% 3150.0	0.045848-0.0027969i
% 3550.0	0.24630-0.24258i
% 4000.0	-0.015962-0.10043i
% 4500.0	-0.053746-0.061193i
% 5000.0	-0.0094247+0.019010i
% ];


%% new values 

folder_data = 'Data';
addpath(folder_data)
data1 = readmatrix('avg0.csv');
data2 = readmatrix('avg6-18.csv');
data3 = readmatrix('avg6-18V.csv');
data4 = readmatrix('avg9-15.csv');
data5 = readmatrix('avg9-15V.csv');
data6 = readmatrix('avg9-27.csv');
data7 = readmatrix('avg9-27V.csv');

%%

Glog1 = 20*log10(abs(data1(:,2)));
Glog2 = 20*log10(abs(data2(:,2)));
Glog3 = 20*log10(abs(data3(:,2)));
Glog4 = 20*log10(abs(data4(:,2)));
Glog5 = 20*log10(abs(data5(:,2)));
Glog6 = 20*log10(abs(data6(:,2)));
Glog7 = 20*log10(abs(data7(:,2)));


%% Plot in semilog
figure(1);
semilogy(data1(:,1), abs(data1(:,2)), 'LineWidth',1.2);
hold on
semilogy(data2(:,1), abs(data2(:,2)), 'LineWidth',1.2);
hold on
semilogy(data3(:,1), abs(data3(:,2)) ,'LineWidth',1.2);
hold on
semilogy(data4(:,1), abs(data4(:,2)),'LineWidth',1.2);
hold on
semilogy(data5(:,1), abs(data5(:,2)),'LineWidth',1.2);
hold on
semilogy(data6(:,1), abs(data6(:,2)),'LineWidth',1.2);
hold on
semilogy(data7(:,1), abs(data7(:,2)),'LineWidth',1.2);
hold off
% title('FR Acoustic Pressure')
lgd = legend('No holes', '6-18', '6-18V', '9-15', '9-15V', '9-27', '9-27V', 'Location','southeast');
set(lgd, 'Box', 'off');
lgd.FontSize = 12;
xlabel('Frequency [Hz]', 'FontSize', 13);
ylabel('Avg Acoustic Pressure [dB]', 'FontSize', 13);
% camroll(90)
% exportgraphics(gcf, 'FRA.pdf'); % use this to save as a cropped pdf

% figure(2);
fh = figure(2);
fh.WindowState = 'maximized';
semilogy(data1(:,1), abs(data1(:,2)), '-.', 'LineWidth',1);
hold on
semilogy(data2(:,1), abs(data2(:,2)), 'LineWidth',1);
hold on
semilogy(data4(:,1), abs(data4(:,2)),'LineWidth',1);
hold on
semilogy(data6(:,1), abs(data6(:,2)),'LineWidth',1);
hold off
% title('FR Acoustic Pressure')
lgd = legend('No holes', '6-18', '9-15', '9-27','Location','southeast');
set(lgd, 'Box', 'off');
lgd.FontSize = 20;
ax = gca;
ax.FontSize = 17; 
xlabel('Frequency [Hz]', 'FontSize', 20);
ylabel('Avg Acoustic Pressure [dB]', 'FontSize', 20);
% xlim([0 1000])
pbaspect([2 1 1])
% exportgraphics(gcf, 'FRAH.pdf'); % use this to save as a cropped pdf



% figure(3);
fh = figure(3);
fh.WindowState = 'maximized';
semilogy(data1(:,1), abs(data1(:,2)), '-.','LineWidth',1);
hold on
semilogy(data3(:,1), abs(data3(:,2)) ,'LineWidth',1);
hold on
semilogy(data5(:,1), abs(data5(:,2)),'LineWidth',1);
hold on
semilogy(data7(:,1), abs(data7(:,2)),'LineWidth',1);
hold off
% title('FR Acoustic Pressure')
lgd = legend('No holes', '6-18V', '9-15V', '9-27V', 'Location','southeast');
set(lgd, 'Box', 'off');
lgd.FontSize = 20;
ax = gca;
ax.FontSize = 17; 
xlabel('Frequency [Hz]', 'FontSize', 20);
ylabel('Avg Acoustic Pressure [dB]', 'FontSize', 20);
xlim([0 1000])
pbaspect([2 1 1])
% exportgraphics(gcf, 'FRAV.pdf'); % use this to save as a cropped pdf




%% Plot in loglog
figure(4);
loglog(data1(:,1), abs(data1(:,2)), 'LineWidth',1.2);
hold on
loglog(data2(:,1), abs(data2(:,2)), 'LineWidth',1.2);
hold on
loglog(data3(:,1), abs(data3(:,2)) ,'LineWidth',1.2);
hold on
loglog(data4(:,1), abs(data4(:,2)),'LineWidth',1.2);
hold on
loglog(data5(:,1), abs(data5(:,2)),'LineWidth',1.2);
hold on
loglog(data6(:,1), abs(data6(:,2)),'LineWidth',1.2);
hold on
loglog(data7(:,1), abs(data7(:,2)),'LineWidth',1.2);
hold off
% title('FR Acoustic Pressure')
lgd = legend('No holes', '6-18', '6-18V', '9-15', '9-15V', '9-27', '9-27V', 'Location','southeast');
set(lgd, 'Box', 'off');
lgd.FontSize = 12;
xlabel('Frequency [Hz]', 'FontSize', 13);
ylabel('Avg Acoustic Pressure [dB]', 'FontSize', 13);
% xlim([0 1000])


% figure(5);
fh = figure(5);
fh.WindowState = 'maximized';
loglog(data1(:,1), abs(data1(:,2)),'--' , 'LineWidth',1.2);
hold on
loglog(data2(:,1), abs(data2(:,2)), 'LineWidth',1.2);
hold on
loglog(data4(:,1), abs(data4(:,2)),'LineWidth',1.2);
hold on
loglog(data6(:,1), abs(data6(:,2)),'LineWidth',1.2);
hold off
% title('FR Acoustic Pressure')
lgd = legend('No holes', '6-18', '9-15', '9-27', 'Location','southeast');
set(lgd, 'Box', 'off');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 13; 
xlabel('Frequency [Hz]', 'FontSize', 14);
ylabel('Avg Acoustic Pressure [dB]', 'FontSize', 14);
xlim([0 1000])
pbaspect([2 1 1])
% exportgraphics(gcf, 'FRAH.pdf'); % use this to save as a cropped pdf

% figure(6);
fh = figure(6);
fh.WindowState = 'maximized';
loglog(data1(:,1), abs(data1(:,2)),'--' , 'LineWidth',1.2);
hold on
loglog(data3(:,1), abs(data3(:,2)) ,'LineWidth',1.2);
hold on
loglog(data5(:,1), abs(data5(:,2)),'LineWidth',1.2);
hold on
loglog(data7(:,1), abs(data7(:,2)),'LineWidth',1.2);
hold off
% title('FR Acoustic Pressure')
lgd = legend('No holes', '6-18V', '9-15V', '9-27V', 'Location','southeast');
set(lgd, 'Box', 'off');
lgd.FontSize = 15;
ax = gca;
ax.FontSize = 13; 
xlabel('Frequency [Hz]', 'FontSize', 14);
ylabel('Avg Acoustic Pressure [dB]', 'FontSize', 14);
xlim([0 1000])
pbaspect([2 1 1])
% exportgraphics(gcf, 'FRAV.pdf'); % use this to save as a cropped pdf
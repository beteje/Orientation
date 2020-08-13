res = 1;
Thresh_circ = 0.9;
W_circ = 7;
W_est = 11;

disp('Single Orientation SNR30');
load('testData_single.mat', 'noisy_data_30')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_30,W_est,W_circ,Thresh_circ,res);
save('Results_single_SNR30.mat','azimuth','elevation','circularity');

disp('Single Orientation SNR20');
load('testData_single.mat', 'noisy_data_20')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_20,W_est,W_circ,Thresh_circ,res);
save('Results_single_SNR20.mat','azimuth','elevation','circularity');

disp('Single Orientation SNR10');
load('testData_single.mat', 'noisy_data_10')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_10,W_est,W_circ,Thresh_circ,res);
save('Results_single_SNR10.mat','azimuth','elevation','circularity');

disp('Double Orientation SNR30');
load('testData_double.mat', 'noisy_data_30')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_30,W_est,W_circ,Thresh_circ,res);
save('Results_double_SNR30.mat','azimuth','elevation','circularity');

disp('Double Orientation SNR20');
load('testData_double.mat', 'noisy_data_20')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_20,W_est,W_circ,Thresh_circ,res);
save('Results_double_SNR20.mat','azimuth','elevation','circularity');

disp('Double Orientation SNR10');
load('testData_double.mat', 'noisy_data_10')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_10,W_est,W_circ,Thresh_circ,res);
save('Results_double_SNR10.mat','azimuth','elevation','circularity');

disp('Random Orientation SNR30')
load('testData_random.mat', 'noisy_data_30')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_30,W_est,W_circ,Thresh_circ,res);
save('Results_random_SNR30.mat','azimuth','elevation','circularity');

disp('Random Orientation SNR20')
load('testData_random.mat', 'noisy_data_20')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_20,W_est,W_circ,Thresh_circ,res);
save('Results_random_SNR20.mat','azimuth','elevation','circularity');

disp('Random Orientation SNR10')
load('testData_random.mat', 'noisy_data_10')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_10,W_est,W_circ,Thresh_circ,res);
save('Results_random_SNR10.mat','azimuth','elevation','circularity');

disp('Range Orientation SNR30')
load('testData_range.mat', 'noisy_data_30')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_30,W_est,W_circ,Thresh_circ,res);
save('Results_range_SNR30.mat','azimuth','elevation','circularity');

disp('Range Orientation SNR20')
load('testData_range.mat', 'noisy_data_20')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_20,W_est,W_circ,Thresh_circ,res);
save('Results_range_SNR20.mat','azimuth','elevation','circularity');

disp('Range Orientation SNR10')
load('testData_range.mat', 'noisy_data_10')
[azimuth,elevation,circularity] = Orient_Est_DS_V2(noisy_data_10,W_est,W_circ,Thresh_circ,res);
save('Results_range_SNR10.mat','azimuth','elevation','circularity');
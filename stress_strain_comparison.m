clear all
close all
clc

% Input - "deform_2.txt" that contains info about the oscillatory
% deformation. Can read multiple such files in the "filename" variable

% Output - 2 text files containing the loss modulus and phase shift.
% Storage modulus can be plotted from these 2 data. 

filename = {'./deform_2.txt'; }; % can call multiple filenames here that 
% can be looped over to perform the analysis in 1 shot. E.g., effect of
% Temperature, frequency, strain amplitude, etc. 

%e.g.:
%filename = {'./300K/deform_2.txt'; './400K/deform_2.txt'; './500K/deform_2.txt'; }; 

freq = 1/ 3 ; % 1/Time period (here, ps^-1). 

N_cycles = 40; % Average data over this many cycles

data_size_repeat = ones(numel(filename),1)*10; % Data in "deform_2.txt" file is written every 10 steps using fix deform command
time_step = 0.002;

for j = 1 : numel(filename)
   
    Data = importdata(filename{j});

    lx = Data.data(:,7);
    ly = Data.data(:,8);
    lz = Data.data(:,9);

    pxy = Data.data(:,10); % in atm - check. 
    pxz = Data.data(:,11);
    pyz = Data.data(:,12);

    xy = Data.data(:,13) - Data.data(1,13); 
    
    shear_strain = xy/lz(1); %% Note - lx = ly = lz in this case!

%     pause
    
    pxy(N_cycles*((1/freq)/time_step)/data_size_repeat(j) + 1 : end) = []; % This removes the data present after "N_cycles" number of cycles. 
    shear_strain(N_cycles*((1/freq)/time_step)/data_size_repeat(j) + 1 : end) = [];

    time_new = (1 : ( numel(pxy)/N_cycles) )*data_size_repeat(j)*time_step;

    indices_1_cycle = (numel(pxy)/N_cycles); % No. of data points in 1 cycle. We are going to average the data across all the shear cycles in a new variable, pxy_new. 

    for i = 1 : indices_1_cycle
        pxy_new(i) = mean(pxy(i : indices_1_cycle : numel(pxy))) * 0.101325; % check conversion
    end

    for i = 1 : indices_1_cycle
        shear_new(i) = mean(shear_strain(i : indices_1_cycle : numel(shear_strain)));
    end

    windowSize = 2;
    plot_variable = pxy_new; 

    % This smoothes out the stress data - may or may not be required depending on how noisy your data is. 
  
%     for i = windowSize/2+1 : windowSize : indices_1_cycle - windowSize/2
% 
%         time_new_new(i) = time_new(i);
%         plot_variable_new(i) = 0.1*mean(plot_variable( (i - windowSize/2) : (i + windowSize/2) )); % in MPa
% 
%     end

% pause

    %%%%%% First get the omega for the strain curve, to be used in fitting
    %%%%%% stress curve. Uses two functions, "curvefit_sin1" and
    %%%%%% "curvefit_sin2" that need to be present in the same folder. 

    [estimates1, model1] = curvefit_sin1(time_new, shear_new);
    y_fit_strain = estimates1(1)*(sin((estimates1(2)*time_new) + estimates1(3)));
    omega = estimates1(2);
    strain_max(j) = estimates1(1);

    %%%%% Now fit a sine curve to the stress curve %%%%

    [estimates2, model2] = curvefit_sin2(time_new, plot_variable, omega);
    y_fit_pressure = estimates2(1)*(sin((omega*time_new) + estimates2(2)));

    stress_max(j) = estimates2(1); % in atm
    phase_diff(j) = estimates2(2)*180/pi; % in degrees

%     if phase_diff(j) < 0
%         phase_diff(j) = phase_diff(j) + 180
%     end

    Modulus(j) = abs(stress_max(j)/strain_max(j)); % in Mpa

    G_loss(j) = Modulus(j)*sin(phase_diff(j)*pi/180);
    
    figure;

    hold on
    box on
    set(gca, 'Linewidth', 2, 'fontsize', 28, 'fontweight', 'bold')

    plot(time_new, plot_variable, '-g', 'LineWidth', 2, 'Markersize', 10);
    plot(time_new, y_fit_pressure, '-b', 'LineWidth', 2, 'Markersize', 10);

    xlabel('Time (ps)', 'FontSize', 28, 'FontWeight', 'Bold');
    ylabel('P_{xy} (MPa)', 'FontSize', 28, 'FontWeight', 'Bold');

    set(gca,'plotboxaspectratio',[1.5 1 1])
    
    str3 = strcat('\delta (degrees) = ', num2str(phase_diff(j)));

    h2 = annotation('textbox', [.2 .4 .1 .1],...
                   'String', str3, 'LineWidth', 2,...
                   'FontSize', 24, 'Fontweight', 'Bold', ...
                   'Tag' , 'annotation');
                      
   clear plot_variable_new plot_variable y_fit_pressure y_fit_strain omega xy shear shear_new pxy pxy_new time time_new;
   
end

figure;

hold on
box on
set(gca, 'Linewidth', 2, 'fontsize', 28, 'fontweight', 'bold')

plot(freq, phase_diff, '-or', 'LineWidth', 2, 'Markersize', 10);

xlabel('Frequency (THz)', 'FontSize', 28, 'FontWeight', 'Bold');
ylabel('Phase difference (degrees)', 'FontSize', 28, 'FontWeight', 'Bold');
set(gca,'plotboxaspectratio',[1.5 1 1])

fid = fopen('freq_phase_diff.dat', 'w');
dlmwrite('freq_phase_diff.dat', [(1:numel(filename))' phase_diff'], 'delimiter', '\t', 'precision', '%.6f')
fclose(fid);

figure;

hold on
box on
set(gca, 'Linewidth', 2, 'fontsize', 28, 'fontweight', 'bold')

plot(freq, G_loss, '-or', 'LineWidth', 2, 'Markersize', 10);

xlabel('Frequency (THz)', 'FontSize', 28, 'FontWeight', 'Bold');
ylabel('Loss modulus (MPa)', 'FontSize', 28, 'FontWeight', 'Bold');
set(gca,'plotboxaspectratio',[1.5 1 1])

fid = fopen('freq_loss_modulus.dat', 'w');
dlmwrite('freq_loss_modulus.dat', [(1:numel(filename))' G_loss'], 'delimiter', '\t', 'precision', '%.6f')
fclose(fid);

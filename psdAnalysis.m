% PDS Analysis code
clear
clc

% Parameters
subjects_normal = {'444','448','455'};
subjects_disease = {'428','447','454'};
subjects_all = cat(2,subjects_normal,subjects_disease);
channel_number = {'MC-R','MC-R','MC-R','MC-R','MC-R','VC-R','VC-L','M-R','MC-R'};

% Plot settings
%--------------------------------------------

% filled regions settings
opacity = 0.3;
color_fill_normal = 0.2*[1,1,1]; % Grey color - adjust coefficient to change darkness
color_fill_disease = 0.2*[1,1,1];

% Default plot settings
color_def = "b";
line_w_def = 1.5;
line_s_def = '-';

% error plot settings
color_error_bar = 'r';                   % color of error bar line
line_w_error = 1.5;                      % thickness of error bar line
marker_face_color = [1,1,1];             % fill color of marker
marker_edge_color = [0, 128, 128] / 255; % color of edge of marker
marker = 'o';                            % marker shape for individual points
marker_size = 7;                         % size of marker for individual points
marker_line_width = 1;                   % marker edge thickness for individual points

% Settings for plots of multiple curves
color_normal= 	"k";       % color for normal mice
line_w_normal = 1;         % line width for normal mice
color_disease = "#0072BD"; % color for diseased mice
line_w_disease = 1;        % line width for diseased mice

% color for brain states
color_wake = "#0072BD";
color_nrem = "#EDB120";
color_rem  = "#A2142F";

% Line style
line_s_wake = '-';
line_s_nrem = '--';
line_s_rem = '-.';

% Font settings
font_name = 'calibri';
font_size = 13;

% set tilting angle for x label
xlabel_angle = 45;

% Set the zoomed-in limits for each brain state
%----------------------------------------------

% Wake zoomed-in limits
xlim_zoom_wake = [0,7];         % X-axis limits for zoomed-in region
ylim_zoom_wake = [1, 1000];     % Y-axis limits for zoomed-in region

% NREM zoomed-in limits
xlim_zoom_nrem = [5,11];        % X-axis limits for zoomed-in region
ylim_zoom_nrem = [1, 100];      % Y-axis limits for zoomed-in region

% REM zoomed-in limits
xlim_zoom_rem = [2,7];          % X-axis limits for zoomed-in region
ylim_zoom_rem = [1000, 300000]; % Y-axis limits for zoomed-in region

% Handle all mice
%---------------------------------------------------------

% Initialise data frame for normal subjects
mice_all = struct();

for i = 1:length(subjects_all)

    % Identify data file corresponding to subject
    psd_wake_file = [subjects_all{i},'_PSD_Wake_data.csv'];
    psd_nrem_file = [subjects_all{i},'_PSD_NREM_data.csv'];
    psd_rem_file = [subjects_all{i},'_PSD_REM_data.csv'];

    mice_all(i).name = subjects_all{i};

    figure
    x_limit_lower = 0;
    x_limit_upper = 50;
    subplot_rows = 4;
    subplot_cols = 2;


    % Wake analysis for particular subject
    %---------------------------------------------------------
    % Extract data for Wake state for particular subject
    mice_all(i).data_wake = readtable(psd_wake_file,VariableNamingRule = 'preserve');

    % Extract data for Wake state for particular subject
    mice_all(i).data_wake = readtable(psd_wake_file,VariableNamingRule = 'preserve');
    mice_all(i).power_wake = table2array(mice_all(i).data_wake(:,2));
    mice_all(i).frequency_wake = table2array(mice_all(i).data_wake(:,1));

    % Wake - normal plot
    subplot(subplot_rows,subplot_cols,1)
    plot(mice_all(i).frequency_wake,mice_all(i).power_wake, Color = color_def, LineStyle=line_s_def, LineWidth = line_w_def)
    grid off; box off
    xlim([x_limit_lower,x_limit_upper])
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Power (μV^2)', FontName = font_name, FontSize = font_size)
    title(['Wake - ',mice_all(i).name], FontName = font_name, FontSize = font_size)

    % Wake - logarithmic plot
    subplot(subplot_rows,subplot_cols,2)
    semilogy(mice_all(i).frequency_wake,mice_all(i).power_wake, Color = color_def, LineStyle=line_s_def, LineWidth=line_w_def)
    grid off; box off
    xlim([x_limit_lower,x_limit_upper])
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size)
    title(['Wake - ',mice_all(i).name], FontName = font_name, FontSize = font_size)

    % NREM analysis for particular subject
    %---------------------------------------------------------
    % Extract data for Wake state for particular subject
    mice_all(i).power_wake = table2array(mice_all(i).data_wake(:,2));

    % Extract data for NREM state for particular subject
    mice_all(i).data_nrem = readtable(psd_nrem_file,VariableNamingRule = 'preserve');
    mice_all(i).power_nrem = table2array(mice_all(i).data_nrem(:,2));
    mice_all(i).frequency_nrem = table2array(mice_all(i).data_nrem(:,1));

    % NREM - normal plot
    subplot(subplot_rows,subplot_cols,3)
    plot(mice_all(i).frequency_nrem,mice_all(i).power_nrem, Color = color_def, LineStyle=line_s_def, LineWidth=line_w_def)
    grid off; box off
    xlim([x_limit_lower,x_limit_upper])
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Power (μV^2)', FontName = font_name, FontSize = font_size)
    title(['NREM - ',mice_all(i).name], FontName = font_name, FontSize = font_size)

    % NREM - logarithmic plot
    subplot(subplot_rows,subplot_cols,4)
    semilogy(mice_all(i).frequency_nrem,mice_all(i).power_nrem, Color = color_def, LineStyle=line_s_def, LineWidth=line_w_def)
    grid off; box off
    xlim([x_limit_lower,x_limit_upper])
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size)
    title(['NREM - ',mice_all(i).name], FontName = font_name, FontSize = font_size)

    % REM analysis for particular subject
    %---------------------------------------------------------
    % Extract data for Wake state for particular subject
    mice_all(i).frequency_wake = table2array(mice_all(i).data_wake(:,1));

    % Extract data for REM state for particular subject
    mice_all(i).data_rem = readtable(psd_rem_file,VariableNamingRule = 'preserve');
    mice_all(i).power_rem = table2array(mice_all(i).data_rem(:,2));
    mice_all(i).frequency_rem = table2array(mice_all(i).data_rem(:,1));

    % REM - normal plot
    subplot(subplot_rows,subplot_cols,5)
    plot(mice_all(i).frequency_rem,mice_all(i).power_rem, Color = color_def, LineStyle=line_s_def, LineWidth=line_w_def)
    grid off; box off
    xlim([x_limit_lower,x_limit_upper])
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Power (μV^2)', FontName = font_name, FontSize = font_size)
    title(['REM - ',mice_all(i).name], FontName = font_name, FontSize = font_size)

    % REM - logarithmic plot
    subplot(subplot_rows,subplot_cols,6)
    semilogy(mice_all(i).frequency_rem,mice_all(i).power_rem, Color = color_def, LineStyle=line_s_def, LineWidth=line_w_def)
    grid off; box off
    xlim([x_limit_lower,x_limit_upper])
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size)
    title(['REM - ',mice_all(i).name], FontName = font_name, FontSize = font_size)

    % Combined Brain States
    %-----------------------------------------------
    % Combined - normal plot
    subplot(subplot_rows,subplot_cols,7); hold on
    plot(mice_all(i).frequency_wake,mice_all(i).power_wake, Color = color_wake, LineStyle=line_s_def, LineWidth=line_w_normal)
    plot(mice_all(i).frequency_nrem,mice_all(i).power_nrem, Color = color_nrem, LineStyle=line_s_def, LineWidth=line_w_normal)
    plot(mice_all(i).frequency_rem,mice_all(i).power_rem, Color = color_rem, LineStyle=line_s_def, LineWidth=line_w_normal)
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Power (μV^2)', FontName = font_name, FontSize = font_size);
    xlim([x_limit_lower,x_limit_upper])
    legend('Wake','NREM','REM', FontName = font_name, FontSize = font_size, Location = 'eastoutside')
    box off

    % Combined - logarithmic plot
    subplot(subplot_rows,subplot_cols,8)
    semilogy(mice_all(i).frequency_wake,mice_all(i).power_wake, Color = color_wake, LineStyle=line_s_def, LineWidth=line_w_normal);hold on
    semilogy(mice_all(i).frequency_nrem,mice_all(i).power_nrem, Color = color_nrem, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
    semilogy(mice_all(i).frequency_rem,mice_all(i).power_rem, Color = color_rem, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
    xlim([x_limit_lower,x_limit_upper])
    legend('Wake','NREM','REM', FontName = font_name, FontSize = font_size, Location = 'eastoutside')
    box off

    sgtitle(subjects_all{i}, FontName = font_name, FontSize = font_size)
    % Close figure
    hold off

    % Generate desired figure for report
    %------------------------------------------------------------
    figure
    % Combined - logarithmic plot
    semilogy(mice_all(i).frequency_wake,mice_all(i).power_wake, Color = color_wake, LineStyle=line_s_def, LineWidth=line_w_normal);hold on
    semilogy(mice_all(i).frequency_nrem,mice_all(i).power_nrem, Color = color_nrem, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
    semilogy(mice_all(i).frequency_rem,mice_all(i).power_rem, Color = color_rem, LineStyle=line_s_def, LineWidth=line_w_normal); hold on

    xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
    ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
    xlim([x_limit_lower,x_limit_upper])

    box off

    legend('Wake','NREM','REM', FontName = font_name, FontSize = font_size, FontWeight = 'bold', Location = 'best')
    title(['Mouse - ',subjects_all{i}], FontWeight = 'bold',FontName = font_name, FontSize = font_size)
end


% Separate Mice into normal and diseased
%---------------------------------------------------------------
normal_mice = mice_all(1:length(subjects_normal));
disease_mice = mice_all(length(subjects_normal)+1:end);

% Create array to store the power_wake values for all normal mice
power_wake_normal = cell2mat({normal_mice(:).power_wake});
power_nrem_normal = cell2mat({normal_mice(:).power_nrem});
power_rem_normal = cell2mat({normal_mice(:).power_rem});

% Calculate mean power for normal mice
avg_normal.power_wake = sum(power_wake_normal,2) / length(normal_mice);
avg_normal.power_nrem = sum(power_nrem_normal,2) / length(normal_mice);
avg_normal.power_rem = sum(power_rem_normal,2) / length(normal_mice);

% Create array to store the power_wake values for all normal mice
power_wake_disease = cell2mat({disease_mice(:).power_wake});
power_nrem_disease = cell2mat({disease_mice(:).power_nrem});
power_rem_disease = cell2mat({disease_mice(:).power_rem});

% Calculate mean power for diseased mice
avg_disease.power_wake = sum(power_wake_disease,2) / length(disease_mice);
avg_disease.power_nrem = sum(power_nrem_disease,2) / length(disease_mice);
avg_disease.power_rem = sum(power_rem_disease,2) / length(disease_mice);

% Average data for all normal subjects and all diseased subjects
%---------------------------------------------------------------
% Since the frequency readings are the same for all mice any one can be
% chosen to represent the average
avg_normal.frequency_wake = normal_mice(1).frequency_wake;
avg_normal.frequency_nrem = normal_mice(1).frequency_nrem;
avg_normal.frequency_rem = normal_mice(1).frequency_rem;

avg_disease.frequency_wake = disease_mice(1).frequency_wake;
avg_disease.frequency_nrem = disease_mice(1).frequency_nrem;
avg_disease.frequency_rem = disease_mice(1).frequency_rem;

% Generate frequency bands information for averaged mice
%--------------------------------------------------------
% Define frequency bands
delta_band = [1, 5];
theta_band = [5, 10];
sigma_band = [10, 16];
gamma_band = [30, 48];

% Calculate errors for averaging plots
%-------------------------------------
% Standard Error is obtained from power readings for all mice of the same
% genotype for each frequency reading

% Normal Mice
%--------------------------------------
avg_wake_normal_all = cell2mat({normal_mice(:).power_wake});
avg_nrem_normal_all = cell2mat({normal_mice(:).power_nrem});
avg_rem_normal_all = cell2mat({normal_mice(:).power_rem});

% Compute standard error for each row
se_poweravg_wake_norm = std(avg_wake_normal_all, 0, 2) / sqrt(size(avg_wake_normal_all, 2));
se_poweravg_nrem_norm = std(avg_nrem_normal_all, 0, 2) / sqrt(size(avg_nrem_normal_all, 2));
se_poweravg_rem_norm = std(avg_rem_normal_all, 0, 2) / sqrt(size(avg_rem_normal_all, 2));

% Diseased Mice
%--------------------------------------
avg_wake_disease_all = cell2mat({disease_mice(:).power_wake});
avg_nrem_disease_all = cell2mat({disease_mice(:).power_nrem});
avg_rem_disease_all = cell2mat({disease_mice(:).power_rem});

% Compute standard error for each row
se_poweravg_wake_disease = std(avg_wake_disease_all, 0, 2) / sqrt(size(avg_wake_disease_all, 2));
se_poweravg_nrem_disease = std(avg_nrem_disease_all, 0, 2) / sqrt(size(avg_nrem_disease_all, 2));
se_poweravg_rem_disease = std(avg_rem_disease_all, 0, 2) / sqrt(size(avg_rem_disease_all, 2));

% Generate plot to show averaging method
%---------------------------------------------------------
% Plot all normal mice for each brain state

% Wake
figure
for i = 1:length(subjects_normal)
    subplot(1,2,1)
    semilogy(normal_mice(i).frequency_wake,normal_mice(i).power_wake, LineStyle = '--'); hold on
end
semilogy(avg_normal.frequency_wake,avg_normal.power_wake,LineStyle = '-', Color = 'k', LineWidth = 1.5)
legend('426-WT' ,'434-WT','444-WT','448-WT','455-WT','Average - WT')
box off
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlim([x_limit_lower,x_limit_upper])
title('WT Averaging', FontWeight = 'bold',FontName = font_name, FontSize = font_size)

for i = 1:length(subjects_disease)
    subplot(1,2,2)
    semilogy(disease_mice(i).frequency_wake,disease_mice(i).power_wake, LineStyle = '--'); hold on
end
semilogy(avg_disease.frequency_wake,avg_disease.power_wake,LineStyle = '-', Color = 'k', LineWidth = 1.5)
legend('428-AS' ,'435-AS','447-AS','454-AS','Average - AS')
box off
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlim([x_limit_lower,x_limit_upper])
title('AS Averaging', FontWeight = 'bold',FontName = font_name, FontSize = font_size)
sgtitle('Wake', FontWeight = 'bold',FontName = font_name, FontSize = font_size)

% NREM
figure
for i = 1:length(subjects_normal)
    subplot(1,2,1)
    semilogy(normal_mice(i).frequency_nrem,normal_mice(i).power_nrem, LineStyle = '--'); hold on
end
semilogy(avg_normal.frequency_nrem,avg_normal.power_nrem,LineStyle = '-', Color = 'k', LineWidth = 1.5)
legend('426-WT' ,'434-WT','444-WT','448-WT','455-WT','Average - WT')
box off
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlim([x_limit_lower,x_limit_upper])
title('WT Averaging', FontWeight = 'bold',FontName = font_name, FontSize = font_size)

for i = 1:length(subjects_disease)
    subplot(1,2,2)
    semilogy(disease_mice(i).frequency_nrem,disease_mice(i).power_nrem, LineStyle = '--'); hold on
end
semilogy(avg_disease.frequency_nrem,avg_disease.power_nrem,LineStyle = '-', Color = 'k', LineWidth = 1.5)
legend('428-AS' ,'435-AS','447-AS','454-AS','Average - AS')
box off
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlim([x_limit_lower,x_limit_upper])
title('AS Averaging', FontWeight = 'bold',FontName = font_name, FontSize = font_size)
sgtitle('NREM', FontWeight = 'bold',FontName = font_name, FontSize = font_size)

% REM
figure
for i = 1:length(subjects_normal)
    subplot(1,2,1)
    semilogy(normal_mice(i).frequency_rem,normal_mice(i).power_rem, LineStyle = '--'); hold on
end
semilogy(avg_normal.frequency_rem,avg_normal.power_rem,LineStyle = '-', Color = 'k', LineWidth = 1.5)
legend('426-WT' ,'434-WT','444-WT','448-WT','455-WT','Average - WT')
box off
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlim([x_limit_lower,x_limit_upper])
title('WT Averaging', FontWeight = 'bold',FontName = font_name, FontSize = font_size)

for i = 1:length(subjects_disease)
    subplot(1,2,2)
    semilogy(disease_mice(i).frequency_rem,disease_mice(i).power_rem, LineStyle = '--'); hold on
end
semilogy(avg_disease.frequency_rem,avg_disease.power_rem,LineStyle = '-', Color = 'k', LineWidth = 1.5)
legend('428-AS' ,'435-AS','447-AS','454-AS','Average - AS')
box off
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlim([x_limit_lower,x_limit_upper])
title('AS Averaging', FontWeight = 'bold',FontName = font_name, FontSize = font_size)
sgtitle('REM', FontWeight = 'bold',FontName = font_name, FontSize = font_size)

% Wake Results
% --------------------------------------------------------------
% Generate plot combining normal and diseased mice - Logarithmic
%---------------------------------------------------------------
figure;

% Generate shading for average power error
%------------------------------------------
error_x = avg_normal.frequency_wake;

norm_error_y1 = avg_normal.power_wake + se_poweravg_wake_norm;
norm_error_y2 = avg_normal.power_wake - se_poweravg_wake_norm;
disease_error_y1 = avg_disease.power_wake + se_poweravg_wake_disease;
disease_error_y2 = avg_disease.power_wake - se_poweravg_wake_disease;

% Generate subplot with brain states of normal and diseased mice
%----------------------------------------------------------------
subplot(2,4,1:2)
semilogy(avg_normal.frequency_wake,avg_normal.power_wake, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_wake,avg_disease.power_wake, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

% Fill the region between errors
fill([error_x; flipud(error_x)], [norm_error_y1; flipud(norm_error_y2)], color_fill_normal, 'FaceAlpha', opacity, 'EdgeColor','none');hold on
fill([error_x; flipud(error_x)], [disease_error_y1; flipud(disease_error_y2)], color_fill_disease, 'FaceAlpha', opacity, 'EdgeColor','none');hold on

% Replot logarithmic curves to be on top
semilogy(avg_normal.frequency_wake,avg_normal.power_wake, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_wake,avg_disease.power_wake, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

grid off; box off
xlim([x_limit_lower,x_limit_upper])
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
legend('Wake WT','Wake AS', AutoUpdate='off',Location='eastoutside', FontName = font_name, FontSize = font_size)

% Generate zoom-in plot of window of interest
%---------------------------------------------
% Create the zoomed-in plot
subplot(2, 4, 3:4);  % 2 rows, 1 column, second plot
semilogy(avg_normal.frequency_wake,avg_normal.power_wake, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_wake,avg_disease.power_wake, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

% Fill the region between errors
fill([error_x; flipud(error_x)], [norm_error_y1; flipud(norm_error_y2)], color_fill_normal, 'FaceAlpha', opacity, 'EdgeColor','none');hold on
fill([error_x; flipud(error_x)], [disease_error_y1; flipud(disease_error_y2)], color_fill_disease, 'FaceAlpha', opacity, 'EdgeColor','none');hold on

% Replot logarithmic curves to be on top
semilogy(avg_normal.frequency_wake,avg_normal.power_wake, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_wake,avg_disease.power_wake, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
legend('Wake WT','Wake AS',Location='eastoutside', FontName = font_name, FontSize = font_size)
grid off; box off

% Set the axis limits for the zoomed-in region
xlim(xlim_zoom_wake);
ylim(ylim_zoom_wake);

% Highlight the zoomed-in region in the original plot
subplot(2, 4, 1:2);  % Switch back to the first subplot
hold on;
plot([xlim_zoom_wake(1), xlim_zoom_wake(1)], ylim_zoom_wake, 'r--');  % Vertical line at the left zoom limit
plot([xlim_zoom_wake(2), xlim_zoom_wake(2)], ylim_zoom_wake, 'r--');  % Vertical line at the right zoom limit
plot(xlim_zoom_wake, [ylim_zoom_wake(1), ylim_zoom_wake(1)], 'r--');  % Horizontal line at the lower zoom limit
plot(xlim_zoom_wake, [ylim_zoom_wake(2), ylim_zoom_wake(2)], 'r--');  % Horizontal line at the upper zoom limit

% Obtain psds (power data) and frequency from overall power spectra on
% selected epochs for normal sibjects
delta_power_band_norm = avg_normal.power_wake((avg_normal.frequency_wake >= delta_band(1)) & (avg_normal.frequency_wake <= delta_band(2)));
theta_power_band_norm = avg_normal.power_wake((avg_normal.frequency_wake >= theta_band(1)) & (avg_normal.frequency_wake <= theta_band(2)));
sigma_power_band_norm = avg_normal.power_wake((avg_normal.frequency_wake >= sigma_band(1)) & (avg_normal.frequency_wake <= sigma_band(2)));
gamma_power_band_norm = avg_normal.power_wake((avg_normal.frequency_wake >= gamma_band(1)) & (avg_normal.frequency_wake <= gamma_band(2)));

% Calculate the mean of obtained power ranges
delta_pmn = mean(delta_power_band_norm); % delta band power mean normal
theta_pmn = mean(theta_power_band_norm);
sigma_pmn = mean(sigma_power_band_norm);
gamma_pmn = mean(gamma_power_band_norm);

% Calculate standard error for normal subjects
se_delta_norm_wake = std(delta_power_band_norm)/sqrt(length(delta_power_band_norm));
se_theta_norm_wake = std(theta_power_band_norm)/sqrt(length(theta_power_band_norm));
se_sigma_norm_wake = std(sigma_power_band_norm)/sqrt(length(sigma_power_band_norm));
se_gamma_norm_wake = std(gamma_power_band_norm)/sqrt(length(gamma_power_band_norm));

% Obtain mean power for diseased subjects
delta_power_band_disease = avg_disease.power_wake((avg_disease.frequency_wake >= delta_band(1)) & (avg_disease.frequency_wake <= delta_band(2)));
theta_power_band_disease = avg_disease.power_wake((avg_disease.frequency_wake >= theta_band(1)) & (avg_disease.frequency_wake <= theta_band(2)));
sigma_power_band_disease = avg_disease.power_wake((avg_disease.frequency_wake >= sigma_band(1)) & (avg_disease.frequency_wake <= sigma_band(2)));
gamma_power_band_disease = avg_disease.power_wake((avg_disease.frequency_wake >= gamma_band(1)) & (avg_disease.frequency_wake <= gamma_band(2)));

% Calculate the mean of obtained power ranges
delta_pmd = mean(delta_power_band_disease); % delta band power mean diseased
theta_pmd = mean(theta_power_band_disease);
sigma_pmd = mean(sigma_power_band_disease);
gamma_pmd = mean(gamma_power_band_disease);

% Calculate standard error for diseased subjects
se_delta_disease_wake = std(delta_power_band_disease)/sqrt(length(delta_power_band_disease));
se_theta_disease_wake = std(theta_power_band_disease)/sqrt(length(theta_power_band_disease));
se_sigma_disease_wake = std(sigma_power_band_disease)/sqrt(length(sigma_power_band_disease));
se_gamma_disease_wake = std(gamma_power_band_disease)/sqrt(length(gamma_power_band_disease));

% Generate bar charts to compare the frequency band power means between
% the genotypes
% =====================================================================
% Delta band
%------------------------------------------------------------
subplot(2,4,5)
delta_bar = bar( [delta_pmn , delta_pmd; 0,0]);hold on % plot bars
delta_bar(1).FaceColor = color_normal;
delta_bar(2).FaceColor = color_disease;
delta_bar(1).EdgeColor = 'none';
delta_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(delta_power_band_norm)), delta_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(delta_power_band_disease)), delta_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],delta_pmn, se_delta_norm_wake, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],delta_pmd, se_delta_disease_wake, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
ylabel('Power (μV^2)', FontName = font_name, FontSize = font_size);
title('Delta', FontName = font_name, FontSize = font_size)

% Theta band
%------------------------------------------------------------
subplot(2,4,6)
theta_bar = bar( [theta_pmn , theta_pmd; 0,0]);hold on
theta_bar(1).FaceColor = color_normal;
theta_bar(2).FaceColor = color_disease;
theta_bar(1).EdgeColor = 'none';
theta_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(theta_power_band_norm)), theta_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(theta_power_band_disease)), theta_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],theta_pmn, se_theta_norm_wake, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],theta_pmd, se_theta_disease_wake, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
title('Theta', FontName = font_name, FontSize = font_size)

% Sigma band
%------------------------------------------------------------
subplot(2,4,7)
sigma_bar = bar( [sigma_pmn , sigma_pmd; 0,0]);hold on
sigma_bar(1).FaceColor = color_normal;
sigma_bar(2).FaceColor = color_disease;
sigma_bar(1).EdgeColor = 'none';
sigma_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(sigma_power_band_norm)), sigma_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(sigma_power_band_disease)), sigma_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],sigma_pmn, se_sigma_norm_wake, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],sigma_pmd, se_sigma_disease_wake, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off;
title('Sigma', FontName = font_name, FontSize = font_size)

% Gamma band
%------------------------------------------------------------
subplot(2,4,8)
gamma_bar = bar( [gamma_pmn , gamma_pmd; 0,0]);hold on
gamma_bar(1).FaceColor = color_normal;
gamma_bar(2).FaceColor = color_disease;
gamma_bar(1).EdgeColor = 'none';
gamma_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(gamma_power_band_norm)), gamma_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(gamma_power_band_disease)), gamma_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],gamma_pmn, se_gamma_norm_wake, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],gamma_pmd, se_gamma_disease_wake, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
title('Gamma', FontName = font_name, FontSize = font_size)
sgtitle("Wake Analysis", FontWeight = 'bold',FontName = font_name, FontSize = font_size)

% NREM Results
% --------------------------------------------------------------
% Generate plot combining normal and diseased mice - Logarithmic
%---------------------------------------------------------------
figure;

% Generate shading for average power error
%------------------------------------------
error_x = avg_normal.frequency_nrem;

norm_error_y1 = avg_normal.power_nrem + se_poweravg_nrem_norm;
norm_error_y2 = avg_normal.power_nrem - se_poweravg_nrem_norm;
disease_error_y1 = avg_disease.power_nrem + se_poweravg_nrem_disease;
disease_error_y2 = avg_disease.power_nrem - se_poweravg_nrem_disease;

% Generate subplot with brain states of normal and diseased mice
subplot(2,4,1:2)
semilogy(avg_normal.frequency_nrem,avg_normal.power_nrem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_nrem,avg_disease.power_nrem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

% Fill the region between errors
fill([error_x; flipud(error_x)], [norm_error_y1; flipud(norm_error_y2)], color_fill_normal, 'FaceAlpha', opacity, 'EdgeColor','none');hold on
fill([error_x; flipud(error_x)], [disease_error_y1; flipud(disease_error_y2)], color_fill_disease, 'FaceAlpha', opacity, 'EdgeColor','none');hold on

% Replot logarithmic curves to be on top
semilogy(avg_normal.frequency_nrem,avg_normal.power_nrem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_nrem,avg_disease.power_nrem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

grid off; box off
xlim([x_limit_lower,x_limit_upper])
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
legend('NREM WT','NREM AS', AutoUpdate='off',Location='eastoutside', FontName = font_name, FontSize = font_size)

% Generate zoom-in plot of window of interest
%---------------------------------------------
% Create the zoomed-in plot
subplot(2, 4, 3:4);  % 2 rows, 1 column, second plot
semilogy(avg_normal.frequency_nrem,avg_normal.power_nrem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_nrem,avg_disease.power_nrem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

% Fill the region between errors
fill([error_x; flipud(error_x)], [norm_error_y1; flipud(norm_error_y2)], color_fill_normal, 'FaceAlpha', opacity, 'EdgeColor','none');hold on
fill([error_x; flipud(error_x)], [disease_error_y1; flipud(disease_error_y2)], color_fill_disease, 'FaceAlpha', opacity, 'EdgeColor','none');hold on

% Replot logarithmic curves to be on top
semilogy(avg_normal.frequency_nrem,avg_normal.power_nrem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_nrem,avg_disease.power_nrem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
legend('NREM WT','NREM AS',Location='eastoutside', FontName = font_name, FontSize = font_size)
grid off; box off

% Set the axis limits for the zoomed-in region
xlim(xlim_zoom_nrem);
ylim(ylim_zoom_nrem);

% Highlight the zoomed-in region in the original plot
subplot(2, 4, 1:2);  % Switch back to the first subplot
hold on;
plot([xlim_zoom_nrem(1), xlim_zoom_nrem(1)], ylim_zoom_nrem, 'r--');  % Vertical line at the left zoom limit
plot([xlim_zoom_nrem(2), xlim_zoom_nrem(2)], ylim_zoom_nrem, 'r--');  % Vertical line at the right zoom limit
plot(xlim_zoom_nrem, [ylim_zoom_nrem(1), ylim_zoom_nrem(1)], 'r--');  % Horizontal line at the lower zoom limit
plot(xlim_zoom_nrem, [ylim_zoom_nrem(2), ylim_zoom_nrem(2)], 'r--');  % Horizontal line at the upper zoom limit

% Obtain psds (power data) and frequency from overall power spectra on
% selected epochs for normal sibjects
delta_power_band_norm = avg_normal.power_nrem((avg_normal.frequency_nrem >= delta_band(1)) & (avg_normal.frequency_nrem <= delta_band(2)));
theta_power_band_norm = avg_normal.power_nrem((avg_normal.frequency_nrem >= theta_band(1)) & (avg_normal.frequency_nrem <= theta_band(2)));
sigma_power_band_norm = avg_normal.power_nrem((avg_normal.frequency_nrem >= sigma_band(1)) & (avg_normal.frequency_nrem <= sigma_band(2)));
gamma_power_band_norm = avg_normal.power_nrem((avg_normal.frequency_nrem >= gamma_band(1)) & (avg_normal.frequency_nrem <= gamma_band(2)));

% Calculate the mean of obtained power ranges
delta_pmn = mean(delta_power_band_norm); % delta band power mean normal
theta_pmn = mean(theta_power_band_norm);
sigma_pmn = mean(sigma_power_band_norm);
gamma_pmn = mean(gamma_power_band_norm);

% Calculate standard error for normal subjects
se_delta_norm_nrem = std(delta_power_band_norm)/sqrt(length(delta_power_band_norm));
se_theta_norm_nrem = std(theta_power_band_norm)/sqrt(length(theta_power_band_norm));
se_sigma_norm_nrem = std(sigma_power_band_norm)/sqrt(length(sigma_power_band_norm));
se_gamma_norm_nrem = std(gamma_power_band_norm)/sqrt(length(gamma_power_band_norm));

% Obtain mean power for diseased subjects
delta_power_band_disease = avg_disease.power_nrem((avg_disease.frequency_nrem >= delta_band(1)) & (avg_disease.frequency_nrem <= delta_band(2)));
theta_power_band_disease = avg_disease.power_nrem((avg_disease.frequency_nrem >= theta_band(1)) & (avg_disease.frequency_nrem <= theta_band(2)));
sigma_power_band_disease = avg_disease.power_nrem((avg_disease.frequency_nrem >= sigma_band(1)) & (avg_disease.frequency_nrem <= sigma_band(2)));
gamma_power_band_disease = avg_disease.power_nrem((avg_disease.frequency_nrem >= gamma_band(1)) & (avg_disease.frequency_nrem <= gamma_band(2)));

% Calculate the mean of obtained power ranges
delta_pmd = mean(delta_power_band_disease); % delta band power mean diseased
theta_pmd = mean(theta_power_band_disease);
sigma_pmd = mean(sigma_power_band_disease);
gamma_pmd = mean(gamma_power_band_disease);

% Calculate standard error for diseased subjects
se_delta_disease_nrem = std(delta_power_band_disease)/sqrt(length(delta_power_band_disease));
se_theta_disease_nrem = std(theta_power_band_disease)/sqrt(length(theta_power_band_disease));
se_sigma_disease_nrem = std(sigma_power_band_disease)/sqrt(length(sigma_power_band_disease));
se_gamma_disease_nrem = std(gamma_power_band_disease)/sqrt(length(gamma_power_band_disease));

% Generate bar charts to compare the frequency band power means between
% the genotypes
% =====================================================================
% Delta band
%------------------------------------------------------------
subplot(2,4,5)
delta_bar = bar( [delta_pmn , delta_pmd; 0,0]);hold on % plot bars
delta_bar(1).FaceColor = color_normal;
delta_bar(2).FaceColor = color_disease;
delta_bar(1).EdgeColor = 'none';
delta_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(delta_power_band_norm)), delta_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(delta_power_band_disease)), delta_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],delta_pmn, se_delta_norm_nrem, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],delta_pmd, se_delta_disease_nrem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
ylabel('Power (μV^2)', FontName = font_name, FontSize = font_size);
title('Delta', FontName = font_name, FontSize = font_size)

% Theta band
%------------------------------------------------------------
subplot(2,4,6)
theta_bar = bar( [theta_pmn , theta_pmd; 0,0]);hold on
theta_bar(1).FaceColor = color_normal;
theta_bar(2).FaceColor = color_disease;
theta_bar(1).EdgeColor = 'none';
theta_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(theta_power_band_norm)), theta_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(theta_power_band_disease)), theta_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],theta_pmn, se_theta_norm_nrem, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],theta_pmd, se_theta_disease_nrem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
title('Theta', FontName = font_name, FontSize = font_size)

% Sigma band
%------------------------------------------------------------
subplot(2,4,7)
sigma_bar = bar( [sigma_pmn , sigma_pmd; 0,0]);hold on
sigma_bar(1).FaceColor = color_normal;
sigma_bar(2).FaceColor = color_disease;
sigma_bar(1).EdgeColor = 'none';
sigma_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(sigma_power_band_norm)), sigma_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(sigma_power_band_disease)), sigma_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],sigma_pmn, se_sigma_norm_nrem, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],sigma_pmd, se_sigma_disease_nrem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off;
title('Sigma', FontName = font_name, FontSize = font_size)

% Gamma band
%------------------------------------------------------------
subplot(2,4,8)
gamma_bar = bar( [gamma_pmn , gamma_pmd; 0,0]);hold on
gamma_bar(1).FaceColor = color_normal;
gamma_bar(2).FaceColor = color_disease;
gamma_bar(1).EdgeColor = 'none';
gamma_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(gamma_power_band_norm)), gamma_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(gamma_power_band_disease)), gamma_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar(0.85,gamma_pmn, se_gamma_norm_nrem, Color = color_error_bar, lineWidth = line_w_error)
errorbar(1.15,gamma_pmd, se_gamma_disease_nrem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
title('Gamma', FontName = font_name, FontSize = font_size)
sgtitle("NREM Analysis", FontWeight = 'bold',FontName = font_name, FontSize = font_size)

% REM Results
% --------------------------------------------------------------
% Generate plot combining normal and diseased mice - Logarithmic
%---------------------------------------------------------------
figure;

% Generate shading for average power error
%------------------------------------------
error_x = avg_normal.frequency_rem;

norm_error_y1 = avg_normal.power_rem + se_poweravg_rem_norm;
norm_error_y2 = avg_normal.power_rem - se_poweravg_rem_norm;
disease_error_y1 = avg_disease.power_rem + se_poweravg_rem_disease;
disease_error_y2 = avg_disease.power_rem - se_poweravg_rem_disease;

% Generate subplot with brain states of normal and diseased mice
subplot(2,4,1:2)
semilogy(avg_normal.frequency_rem,avg_normal.power_rem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_rem,avg_disease.power_rem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

% Fill the region between errors
fill([error_x; flipud(error_x)], [norm_error_y1; flipud(norm_error_y2)], color_fill_normal, 'FaceAlpha', opacity, 'EdgeColor','none');hold on
fill([error_x; flipud(error_x)], [disease_error_y1; flipud(disease_error_y2)], color_fill_disease, 'FaceAlpha', opacity, 'EdgeColor','none');hold on

% Replot logarithmic curves to be on top
semilogy(avg_normal.frequency_rem,avg_normal.power_rem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_rem,avg_disease.power_rem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

grid off; box off
xlim([x_limit_lower,x_limit_upper])
ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
legend('REM WT','REM AS', AutoUpdate='off',Location='eastoutside', FontName = font_name, FontSize = font_size)

% Generate zoom-in plot of window of interest
%---------------------------------------------

% Create the zoomed-in plot
subplot(2, 4, 3:4);  % 2 rows, 1 column, second plot
semilogy(avg_normal.frequency_rem,avg_normal.power_rem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_rem,avg_disease.power_rem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

% Fill the region between errors
fill([error_x; flipud(error_x)], [norm_error_y1; flipud(norm_error_y2)], color_fill_normal, 'FaceAlpha', opacity, 'EdgeColor','none');hold on
fill([error_x; flipud(error_x)], [disease_error_y1; flipud(disease_error_y2)], color_fill_disease, 'FaceAlpha', opacity, 'EdgeColor','none');hold on

% Replot logarithmic curves to be on top
semilogy(avg_normal.frequency_rem,avg_normal.power_rem, Color = color_normal, LineStyle=line_s_def, LineWidth=line_w_normal); hold on
semilogy(avg_disease.frequency_rem,avg_disease.power_rem, Color = color_disease, LineStyle=line_s_def, LineWidth=line_w_disease); hold on

ylabel('Log Power (μV^2)', FontName = font_name, FontSize = font_size);
xlabel('Frequency (Hz)', FontName = font_name, FontSize = font_size)
legend('REM WT','REM AS',Location='eastoutside', FontName = font_name, FontSize = font_size)
grid off; box off

% Set the axis limits for the zoomed-in region
xlim(xlim_zoom_rem);
ylim(ylim_zoom_rem);

% Highlight the zoomed-in region in the original plot
subplot(2, 4, 1:2);  % Switch back to the first subplot
hold on;
plot([xlim_zoom_rem(1), xlim_zoom_rem(1)], ylim_zoom_rem, 'r--');  % Vertical line at the left zoom limit
plot([xlim_zoom_rem(2), xlim_zoom_rem(2)], ylim_zoom_rem, 'r--');  % Vertical line at the right zoom limit
plot(xlim_zoom_rem, [ylim_zoom_rem(1), ylim_zoom_rem(1)], 'r--');  % Horizontal line at the lower zoom limit
plot(xlim_zoom_rem, [ylim_zoom_rem(2), ylim_zoom_rem(2)], 'r--');  % Horizontal line at the upper zoom limit

% Obtain psds (power data) and frequency from overall power spectra on
% selected epochs for normal sibjects
delta_power_band_norm = avg_normal.power_rem((avg_normal.frequency_rem >= delta_band(1)) & (avg_normal.frequency_rem <= delta_band(2)));
theta_power_band_norm = avg_normal.power_rem((avg_normal.frequency_rem >= theta_band(1)) & (avg_normal.frequency_rem <= theta_band(2)));
sigma_power_band_norm = avg_normal.power_rem((avg_normal.frequency_rem >= sigma_band(1)) & (avg_normal.frequency_rem <= sigma_band(2)));
gamma_power_band_norm = avg_normal.power_rem((avg_normal.frequency_rem >= gamma_band(1)) & (avg_normal.frequency_rem <= gamma_band(2)));

% Calculate the mean of obtained power ranges
delta_pmn = mean(delta_power_band_norm); % delta band power mean normal
theta_pmn = mean(theta_power_band_norm);
sigma_pmn = mean(sigma_power_band_norm);
gamma_pmn = mean(gamma_power_band_norm);

% Calculate standard error for normal subjects
se_delta_norm_rem = std(delta_power_band_norm)/sqrt(length(delta_power_band_norm));
se_theta_norm_rem = std(theta_power_band_norm)/sqrt(length(theta_power_band_norm));
se_sigma_norm_rem = std(sigma_power_band_norm)/sqrt(length(sigma_power_band_norm));
se_gamma_norm_rem = std(gamma_power_band_norm)/sqrt(length(gamma_power_band_norm));

% Obtain mean power for diseased subjects
delta_power_band_disease = avg_disease.power_rem((avg_disease.frequency_rem >= delta_band(1)) & (avg_disease.frequency_rem <= delta_band(2)));
theta_power_band_disease = avg_disease.power_rem((avg_disease.frequency_rem >= theta_band(1)) & (avg_disease.frequency_rem <= theta_band(2)));
sigma_power_band_disease = avg_disease.power_rem((avg_disease.frequency_rem >= sigma_band(1)) & (avg_disease.frequency_rem <= sigma_band(2)));
gamma_power_band_disease = avg_disease.power_rem((avg_disease.frequency_rem >= gamma_band(1)) & (avg_disease.frequency_rem <= gamma_band(2)));

% Calculate the mean of obtained power ranges
delta_pmd = mean(delta_power_band_disease); % delta band power mean diseased
theta_pmd = mean(theta_power_band_disease);
sigma_pmd = mean(sigma_power_band_disease);
gamma_pmd = mean(gamma_power_band_disease);

% Calculate standard error for diseased subjects
se_delta_disease_rem = std(delta_power_band_disease)/sqrt(length(delta_power_band_disease));
se_theta_disease_rem = std(theta_power_band_disease)/sqrt(length(theta_power_band_disease));
se_sigma_disease_rem = std(sigma_power_band_disease)/sqrt(length(sigma_power_band_disease));
se_gamma_disease_rem = std(gamma_power_band_disease)/sqrt(length(gamma_power_band_disease));

% Generate bar charts to compare the frequency band power means between
% the genotypes
% =====================================================================
% Delta band
%------------------------------------------------------------
subplot(2,4,5)
delta_bar = bar( [delta_pmn , delta_pmd; 0,0]);hold on % plot bars
delta_bar(1).FaceColor = color_normal;
delta_bar(2).FaceColor = color_disease;
delta_bar(1).EdgeColor = 'none';
delta_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(delta_power_band_norm)), delta_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(delta_power_band_disease)), delta_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],delta_pmn, se_delta_norm_rem, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],delta_pmd, se_delta_disease_rem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
ylabel('Power (μV^2)', FontName = font_name, FontSize = font_size);
title('Delta', FontName = font_name, FontSize = font_size)

% Theta band
%------------------------------------------------------------
subplot(2,4,6)
theta_bar = bar( [theta_pmn , theta_pmd; 0,0]);hold on
theta_bar(1).FaceColor = color_normal;
theta_bar(2).FaceColor = color_disease;
theta_bar(1).EdgeColor = 'none';
theta_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(theta_power_band_norm)), theta_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(theta_power_band_disease)), theta_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],theta_pmn, se_theta_norm_rem, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],theta_pmd, se_theta_disease_rem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
title('Theta', FontName = font_name, FontSize = font_size)

% Sigma band
%------------------------------------------------------------
subplot(2,4,7)
sigma_bar = bar( [sigma_pmn , sigma_pmd; 0,0]);hold on
sigma_bar(1).FaceColor = color_normal;
sigma_bar(2).FaceColor = color_disease;
sigma_bar(1).EdgeColor = 'none';
sigma_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(sigma_power_band_norm)), sigma_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(sigma_power_band_disease)), sigma_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],sigma_pmn, se_sigma_norm_rem, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],sigma_pmd, se_sigma_disease_rem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off;
title('Sigma', FontName = font_name, FontSize = font_size)

% Gamma band
%------------------------------------------------------------
subplot(2,4,8)
gamma_bar = bar( [gamma_pmn , gamma_pmd; 0,0]);hold on
gamma_bar(1).FaceColor = color_normal;
gamma_bar(2).FaceColor = color_disease;
gamma_bar(1).EdgeColor = 'none';
gamma_bar(2).EdgeColor = 'none';

% Plot individual readings
plot(0.85 * ones(length(gamma_power_band_norm)), gamma_power_band_norm,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
plot(1.15 * ones(length(gamma_power_band_disease)), gamma_power_band_disease,Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

% Plot standard error
errorbar([0.85,0.85],gamma_pmn, se_gamma_norm_rem, Color = color_error_bar, lineWidth = line_w_error)
errorbar([1.15,1.15],gamma_pmd, se_gamma_disease_rem, Color = color_error_bar, lineWidth = line_w_error)

xlim([0.5,1.5])
xticks([0.85,1.15])
xticklabels({'WT','HET'})
set(gca, FontName = font_name, FontSize = font_size)
xtickangle(xlabel_angle)

box off
title('Gamma', FontName = font_name, FontSize = font_size)
sgtitle("REM Analysis", FontWeight = 'bold',FontName = font_name, FontSize = font_size)

% Collect Standard Errors for all frequency bands
SE_all = [se_delta_norm_wake,se_delta_disease_wake,se_theta_norm_wake,se_theta_disease_wake,se_sigma_norm_wake,se_sigma_disease_wake,se_gamma_norm_wake,se_gamma_disease_wake;
    se_delta_norm_nrem,se_delta_disease_nrem,se_theta_norm_nrem,se_theta_disease_nrem,se_sigma_norm_nrem,se_sigma_disease_nrem,se_gamma_norm_nrem,se_gamma_disease_nrem;
    se_delta_norm_rem,se_delta_disease_rem,se_theta_norm_rem,se_theta_disease_rem,se_sigma_norm_rem,se_sigma_disease_rem,se_gamma_norm_rem,se_gamma_disease_rem;]';

row_names = {'Delta - WT','Delta - AT','Theta - WT','Theta - AT','Sigma - WT','Sigma - AT','Gamma - WT','Gamma - AT'};
column_names = {'Wake','NREM','REM'};
Table_SE = array2table(SE_all,'VariableNames',column_names,'RowNames',row_names);
disp(Table_SE)

% Create a summary table for all the data

% Initialize the cell array
psd_summary_temp = cell(length(subjects_all), 1);

% Create a summary table for all the data
for i = 1:length(subjects_all)
    if i <= length(subjects_normal)
        genotype = 'WT';
    else
        genotype = 'AS';
    end

    psd_summary_name    = repmat({subjects_all{i}}, length(mice_all(i).power_wake), 1);
    psd_summary_channel = repmat({channel_number{i}}, length(mice_all(i).power_wake), 1);

    % Handle power and frequency data
    psd_summary_power_wake     = mice_all(i).power_wake;
    psd_summary_frequency_wake = mice_all(i).frequency_wake;
    psd_summary_power_nrem     = mice_all(i).power_nrem;
    psd_summary_frequency_nrem = mice_all(i).frequency_nrem;
    psd_summary_power_rem      = mice_all(i).power_rem;
    psd_summary_frequency_rem  = mice_all(i).frequency_rem;

    % Handle frequency bands
    psd_summary_freq_band = cell(length(mice_all(i).frequency_wake),1);
    psd_summary_freq_band((mice_all(i).frequency_wake < delta_band(1)),1) = {'other'};
    psd_summary_freq_band((mice_all(i).frequency_wake >= delta_band(1)) & (mice_all(i).frequency_wake < delta_band(2)),1) = {'delta'};
    psd_summary_freq_band((mice_all(i).frequency_wake >= theta_band(1)) & (mice_all(i).frequency_wake < theta_band(2)),1) = {'theta'};
    psd_summary_freq_band((mice_all(i).frequency_wake >= sigma_band(1)) & (mice_all(i).frequency_wake < sigma_band(2)),1) = {'sigma'};
    psd_summary_freq_band((mice_all(i).frequency_wake > sigma_band(2)) & (mice_all(i).frequency_wake <= gamma_band(1)),1) = {'other'};
    psd_summary_freq_band((mice_all(i).frequency_wake >= gamma_band(1)) & (mice_all(i).frequency_wake <= gamma_band(2)),1) = {'gamma'};
    psd_summary_freq_band((mice_all(i).frequency_wake > gamma_band(2)),1) = {'other'};

    psd_summary_genotype = repmat({genotype}, length(mice_all(i).power_wake), 1);

    % Combine in summary monomer
    psd_summary_temp{i, 1} = {psd_summary_name, psd_summary_channel, psd_summary_power_wake,psd_summary_frequency_wake,psd_summary_freq_band,psd_summary_power_nrem,psd_summary_frequency_nrem,psd_summary_freq_band, psd_summary_power_rem,psd_summary_frequency_rem,psd_summary_freq_band, psd_summary_genotype};
end

merged_summary = vertcat(psd_summary_temp{:});

summary_animal_id = vertcat(merged_summary{:,1});
summary_channel_num = vertcat(merged_summary{:,2});
summary_power_wake = vertcat(merged_summary{:,3});
summary_frequency_wake = vertcat(merged_summary{:,4});
summary_power_nrem = vertcat(merged_summary{:,6});
summary_frequency_nrem = vertcat(merged_summary{:,7});
summary_power_rem = vertcat(merged_summary{:,9});
summary_frequency_rem = vertcat(merged_summary{:,10});
summary_genotype = vertcat(merged_summary{:,12});
summary_frequency_band = vertcat(merged_summary{:,5});

psd_summary = {summary_animal_id,summary_channel_num,summary_power_wake,summary_frequency_wake,summary_frequency_band, summary_power_nrem,summary_frequency_nrem,summary_frequency_band, summary_power_rem,summary_frequency_rem,summary_frequency_band, summary_genotype};
T_summary = table(psd_summary{1},psd_summary{2},psd_summary{3},psd_summary{4},psd_summary{5},psd_summary{6},psd_summary{7},psd_summary{8},psd_summary{9},psd_summary{10},psd_summary{11},psd_summary{12},'VariableNames',{'Animal ID','Channel','Power - Wake','Frequency - Wake','Frequency Band - Wake','Power - NREM','Frequency - NREM','Frequency Band - NREM','Power - REM','Frequency - REM','Frequency Band - REM','Genotype'});

% Display the table
disp(T_summary)

% Export the table to a CSV file
writetable(T_summary, 'PSD_summary_filtered.csv');

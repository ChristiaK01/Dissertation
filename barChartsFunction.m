% Generates barcharts for any data
clear
clc


% Plot settings
%--------------------------------------------
plot_set = struct();


% Default plot settings
plot_set.color_def = "b";
plot_set.line_w_def = 1.5;
plot_set.line_s_def = '-';

% error plot settings
plot_set.color_error_bar = 'r';                   % color of error bar line
plot_set.line_w_error = 1.5;                      % thickness of error bar line
plot_set.marker_face_color = [1,1,1];             % fill color of marker
plot_set.marker_edge_color = [0, 128, 128] / 255; % color of edge of marker
plot_set.marker = 'o';                            % marker shape for individual points
plot_set.marker_size = 7;                         % size of marker for individual points
plot_set.marker_line_width = 1;                   % marker edge thickness for individual points

% Settings for plots of multiple curves
plot_set.color_wt= 	"k";       % color for normal mice
plot_set.line_w_wt = 1;         % line width for normal mice
plot_set.color_het = "#0072BD"; % color for diseased mice
plot_set.line_w_het = 1;        % line width for diseased mice

% color for brain states
plot_set.color_wake = "#0072BD";
plot_set.color_nrem = "#EDB120";
plot_set.color_rem  = "#A2142F";

% Line style
plot_set.line_s_wake = '-';
plot_set.line_s_nrem = '--';
plot_set.line_s_rem = '-.';

% Font settings
plot_set.font_name = 'calibri';
plot_set.font_size = 13;
plot_set.font_size_title = 15;

% set tilting angle for x label
plot_set.xlabel_angle = 45;

% y labels for each barchart
plot_set.ylabels = {'Power (μV^2)'};

% Organise Data
% ==================================================================
power.data_all = readtable("data_all_power.xlsx",'VariableNamingRule','preserve');
stats.data_all = readtable("stats_results.xlsx",'ReadRowNames',true,'VariableNamingRule','preserve');
stats.power_avg = stats.data_all.('Avg power');
stats.power_sem = stats.data_all.('SEM');
stats.p_value = stats.data_all.('p');

% WT data
%-------------------------------------------------------------------
% delta
power.wt.delta.wake = power.data_all.('wake-delta-wt');
power.wt.delta.nrem = power.data_all.('nrem-delta-wt');
power.wt.delta.rem = power.data_all.('rem-delta-wt');

stats.wt.delta.wake.power_avg = stats.power_avg(1);
stats.wt.delta.nrem.power_avg = stats.power_avg(5);
stats.wt.delta.rem.power_avg  = stats.power_avg(9);

stats.wt.delta.wake.sem = stats.power_sem(1);
stats.wt.delta.nrem.sem = stats.power_sem(5);
stats.wt.delta.rem.sem  = stats.power_sem(9);

% theta
power.wt.theta.wake = power.data_all.('wake-theta-wt');
power.wt.theta.nrem = power.data_all.('nrem-theta-wt');
power.wt.theta.rem = power.data_all.('rem-theta-wt');

stats.wt.theta.wake.power_avg = stats.power_avg(2);
stats.wt.theta.nrem.power_avg = stats.power_avg(6);
stats.wt.theta.rem.power_avg  = stats.power_avg(10);

stats.wt.theta.wake.sem = stats.power_sem(2);
stats.wt.theta.nrem.sem = stats.power_sem(6);
stats.wt.theta.rem.sem  = stats.power_sem(10);

% sigma
power.wt.sigma.wake = power.data_all.('wake-sigma-wt');
power.wt.sigma.nrem = power.data_all.('nrem-sigma-wt');
power.wt.sigma.rem = power.data_all.('rem-sigma-wt');

stats.wt.sigma.wake.power_avg = stats.power_avg(3);
stats.wt.sigma.nrem.power_avg = stats.power_avg(7);
stats.wt.sigma.rem.power_avg  = stats.power_avg(11);

stats.wt.sigma.wake.sem = stats.power_sem(3);
stats.wt.sigma.nrem.sem = stats.power_sem(7);
stats.wt.sigma.rem.sem  = stats.power_sem(11);

% gamma
power.wt.gamma.wake = power.data_all.('wake-gamma-wt');
power.wt.gamma.nrem = power.data_all.('nrem-gamma-wt');
power.wt.gamma.rem = power.data_all.('rem-gamma-wt');

stats.wt.gamma.wake.power_avg = stats.power_avg(4);
stats.wt.gamma.nrem.power_avg = stats.power_avg(8);
stats.wt.gamma.rem.power_avg  = stats.power_avg(12);

stats.wt.gamma.wake.sem = stats.power_sem(4);
stats.wt.gamma.nrem.sem = stats.power_sem(8);
stats.wt.gamma.rem.sem  = stats.power_sem(12);

% HET data
%-------------------------------------------------------------------
% delta
power.het.delta.wake = power.data_all.('wake-delta-het');
power.het.delta.nrem = power.data_all.('nrem-delta-het');
power.het.delta.rem = power.data_all.('rem-delta-het');

stats.het.delta.wake.power_avg = stats.power_avg(13);
stats.het.delta.nrem.power_avg = stats.power_avg(17);
stats.het.delta.rem.power_avg  = stats.power_avg(21);

stats.het.delta.wake.sem = stats.power_sem(13);
stats.het.delta.nrem.sem = stats.power_sem(17);
stats.het.delta.rem.sem  = stats.power_sem(21);

% theta
power.het.theta.wake = power.data_all.('wake-theta-het');
power.het.theta.nrem = power.data_all.('nrem-theta-het');
power.het.theta.rem = power.data_all.('rem-theta-het');

stats.het.theta.wake.power_avg = stats.power_avg(14);
stats.het.theta.nrem.power_avg = stats.power_avg(18);
stats.het.theta.rem.power_avg  = stats.power_avg(22);

stats.het.theta.wake.sem = stats.power_sem(14);
stats.het.theta.nrem.sem = stats.power_sem(18);
stats.het.theta.rem.sem  = stats.power_sem(22);

% sigma
power.het.sigma.wake = power.data_all.('wake-sigma-het');
power.het.sigma.nrem = power.data_all.('nrem-sigma-het');
power.het.sigma.rem = power.data_all.('rem-sigma-het');

stats.het.sigma.wake.power_avg = stats.power_avg(15);
stats.het.sigma.nrem.power_avg = stats.power_avg(19);
stats.het.sigma.rem.power_avg  = stats.power_avg(23);

stats.het.sigma.wake.sem = stats.power_sem(15);
stats.het.sigma.nrem.sem = stats.power_sem(19);
stats.het.sigma.rem.sem  = stats.power_sem(23);

% gamma
power.het.gamma.wake = power.data_all.('wake-gamma-het');
power.het.gamma.nrem = power.data_all.('nrem-gamma-het');
power.het.gamma.rem = power.data_all.('rem-gamma-het');

stats.het.gamma.wake.power_avg = stats.power_avg(16);
stats.het.gamma.nrem.power_avg = stats.power_avg(20);
stats.het.gamma.rem.power_avg  = stats.power_avg(24);

stats.het.gamma.wake.sem = stats.power_sem(16);
stats.het.gamma.nrem.sem = stats.power_sem(20);
stats.het.gamma.rem.sem  = stats.power_sem(24);

% Generate Wake Bar Charts
% =========================================================================

% Collect data for function input
%---------------------------------
% collect the values for the bars
bar_values = [stats.wt.delta.wake.power_avg,stats.het.delta.wake.power_avg;
    stats.wt.theta.wake.power_avg,stats.het.theta.wake.power_avg;
    stats.wt.sigma.wake.power_avg,stats.het.sigma.wake.power_avg;
    stats.wt.gamma.wake.power_avg,stats.het.gamma.wake.power_avg];

% collect the individual points to be displayed over the bars
individual_points_1 = {power.wt.delta.wake, power.wt.theta.wake, power.wt.sigma.wake, power.wt.gamma.wake};
individual_points_2 = {power.het.delta.wake, power.het.theta.wake, power.het.sigma.wake, power.het.gamma.wake};

% collect the standard errors to use for errorbars
sem_1 = [stats.wt.delta.wake.sem, stats.wt.theta.wake.sem, stats.wt.sigma.wake.sem, stats.wt.gamma.wake.sem];
sem_2 = [stats.het.delta.wake.sem, stats.het.theta.wake.sem, stats.het.sigma.wake.sem, stats.het.gamma.wake.sem];

% Generate bar charts
% --------------------------------
figure
wake_barchart = barcharts(bar_values, individual_points_1, individual_points_2, sem_1, sem_2, plot_set);
sgtitle('Wake', FontWeight = 'bold', FontName = plot_set.font_name, FontSize = plot_set.font_size_title)
hold off

% Generate NREM Bar Charts
% =========================================================================

% Collect data for function input
%---------------------------------
% collect the values for the bars
bar_values = [stats.wt.delta.nrem.power_avg,stats.het.delta.nrem.power_avg;
    stats.wt.theta.nrem.power_avg,stats.het.theta.nrem.power_avg;
    stats.wt.sigma.nrem.power_avg,stats.het.sigma.nrem.power_avg;
    stats.wt.gamma.nrem.power_avg,stats.het.gamma.nrem.power_avg];

% collect the individual points to be displayed over the bars
individual_points_1 = {power.wt.delta.nrem, power.wt.theta.nrem, power.wt.sigma.nrem, power.wt.gamma.nrem};
individual_points_2 = {power.het.delta.nrem, power.het.theta.nrem, power.het.sigma.nrem, power.het.gamma.nrem};

% collect the standard errors to use for errorbars
sem_1 = [stats.wt.delta.nrem.sem, stats.wt.theta.nrem.sem, stats.wt.sigma.nrem.sem, stats.wt.gamma.nrem.sem];
sem_2 = [stats.het.delta.nrem.sem, stats.het.theta.nrem.sem, stats.het.sigma.nrem.sem, stats.het.gamma.nrem.sem];

% Generate bar charts
% --------------------------------
figure
nrem_barchart = barcharts(bar_values, individual_points_1, individual_points_2, sem_1, sem_2, plot_set);
sgtitle('NREM', FontWeight = 'bold', FontName = plot_set.font_name, FontSize = plot_set.font_size_title)
hold off

% Generate REM Bar Charts
% =========================================================================

% Collect data for function input
%---------------------------------
% collect the values for the bars
bar_values = [stats.wt.delta.rem.power_avg,stats.het.delta.rem.power_avg;
    stats.wt.theta.rem.power_avg,stats.het.theta.rem.power_avg;
    stats.wt.sigma.rem.power_avg,stats.het.sigma.rem.power_avg;
    stats.wt.gamma.rem.power_avg,stats.het.gamma.rem.power_avg];

% collect the individual points to be displayed over the bars
individual_points_1 = {power.wt.delta.rem, power.wt.theta.rem, power.wt.sigma.rem, power.wt.gamma.rem};
individual_points_2 = {power.het.delta.rem, power.het.theta.rem, power.het.sigma.rem, power.het.gamma.rem};

% collect the standard errors to use for errorbars
sem_1 = [stats.wt.delta.rem.sem, stats.wt.theta.rem.sem, stats.wt.sigma.rem.sem, stats.wt.gamma.rem.sem];
sem_2 = [stats.het.delta.rem.sem, stats.het.theta.rem.sem, stats.het.sigma.rem.sem, stats.het.gamma.rem.sem];

% Generate bar charts
% --------------------------------
figure
rem_barchart = barcharts(bar_values, individual_points_1, individual_points_2, sem_1, sem_2, plot_set);
sgtitle('REM', FontWeight = 'bold', FontName = plot_set.font_name, FontSize = plot_set.font_size_title)
hold off

function [barchart] = barcharts(bar_values, individual_points_1, individual_points_2, sem_1, sem_2, plot_set)

subplot_title = {'Delta','Theta','Sigma','Gamma'};
for i = 1 : 4
    subplot(1,4,i); hold on

    % Plot bars
    barchart = bar([bar_values(i,:);0,0]);
    % Set plot settings
    barchart(1).FaceColor = plot_set.color_wt;
    barchart(2).FaceColor = plot_set.color_het;
    barchart(1).EdgeColor = 'none';
    barchart(2).EdgeColor = 'none';

    % Plot individual readings
    plot(0.85 * ones(length(individual_points_1{i})), individual_points_1{i}, Marker = plot_set.marker, MarkerFaceColor = plot_set.marker_face_color,MarkerEdgeColor = plot_set.marker_edge_color,LineStyle = 'none', MarkerSize = plot_set.marker_size, LineWidth = plot_set.marker_line_width)
    plot(1.15 * ones(length(individual_points_2{i})), individual_points_2{i}, Marker = plot_set.marker, MarkerFaceColor = plot_set.marker_face_color,MarkerEdgeColor = plot_set.marker_edge_color,LineStyle = 'none', MarkerSize = plot_set.marker_size, LineWidth = plot_set.marker_line_width)

    % Plot standard error
    errorbar(0.85, bar_values(i,1), sem_1, Color = plot_set.color_error_bar, lineWidth = plot_set.line_w_error)
    errorbar(1.15, bar_values(i,2), sem_2, Color = plot_set.color_error_bar, lineWidth = plot_set.line_w_error)

    xlim([0.5,1.5])
    box off
    ylabel(plot_set.ylabels, FontName = plot_set.font_name, FontSize = plot_set.font_size);
    xticks([0.85,1.15])
    xticklabels({'WT','HET'})
    xtickangle(plot_set.xlabel_angle)
    set(gca, FontName = plot_set.font_name, FontSize = plot_set.font_size)
    title(subplot_title{i}, FontName = plot_set.font_name, FontSize = plot_set.font_size)
end

end

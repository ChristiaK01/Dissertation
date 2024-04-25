% Bar Charts Plotting - BOUT
% ====================================================================
clear; clc

files = {"percentage.csv","boutnumber.csv","boutduration.csv"};

% Plot settings
%--------------------------------------------

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
color_wt= 	"k";       % color for normal mice
line_w_wt = 1;         % line width for normal mice
color_het = "#0072BD"; % color for diseased mice
line_w_het = 1;        % line width for diseased mice

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
font_size_title = 15;

% set tilting angle for x label
xlabel_angle = 45;

% y labels for each barchart
ylabels = {'Time in state (%)','Number of Bouts','Average Bout Duration (s)'};

% Obtain data from each file
sets_all = struct();
test_het = struct();
test_wt = struct();

for i = 1:length(files)

    sets_all(i).data = table2cell(readtable(files{i},'VariableNamingRule','preserve'));

    % Separate genotypes
    het_idx = find(strcmp(sets_all(i).data(:,2),'HET')); % find indecies of HET genotype
    wt_idx  = find(strcmp(sets_all(i).data(:,2),'WT'));  % find indecies of WT genotype
    test_het(i).data = sets_all(i).data(het_idx,:);         % isolate HET subjects
    test_wt(i).data  = sets_all(i).data(wt_idx,:);          % isolate WT subjects

    % Obtain each brain state for HET mice
    test_het(i).wake_all = cell2mat(test_het(i).data(:,4));
    test_het(i).nrem_all = cell2mat(test_het(i).data(:,5));
    test_het(i).rem_all  = cell2mat(test_het(i).data(:,6));
    test_het(i).subjects = test_het(i).data(:,1);

    % Obtain each brain state for WT mice
    test_wt(i).wake_all = cell2mat(test_wt(i).data(:,4));
    test_wt(i).nrem_all = cell2mat(test_wt(i).data(:,5));
    test_wt(i).rem_all  = cell2mat(test_wt(i).data(:,6));
    test_wt(i).subjects = test_wt(i).data(:,1);

    % Calculate mean and SEM for HET for each brain state
    test_het(i).wake_avg = mean(test_het(i).wake_all);
    test_het(i).wake_sem = std(test_het(i).wake_all) / sqrt(length(test_het(i).wake_all));
    test_het(i).nrem_avg = mean(test_het(i).nrem_all);
    test_het(i).nrem_sem = std(test_het(i).nrem_all) / sqrt(length(test_het(i).nrem_all));
    test_het(i).rem_avg  = mean(test_het(i).rem_all);
    test_het(i).rem_sem = std(test_het(i).rem_all) / sqrt(length(test_het(i).rem_all));

    % Calculate mean and SEM for WT for each brain state
    test_wt(i).wake_avg = mean(test_wt(i).wake_all);
    test_wt(i).wake_sem = std(test_wt(i).wake_all) / sqrt(length(test_wt(i).wake_all));
    test_wt(i).nrem_avg = mean(test_wt(i).nrem_all);
    test_wt(i).nrem_sem = std(test_wt(i).nrem_all) / sqrt(length(test_wt(i).nrem_all));
    test_wt(i).rem_avg  = mean(test_wt(i).rem_all);
    test_wt(i).rem_sem = std(test_wt(i).rem_all) / sqrt(length(test_wt(i).rem_all));
end


% Generate Bar charts
% ================================================================================

%---------------------------------------------------------------------------------
% Wake
%---------------------------------------------------------------------------------
figure
for i = 1:length(files)
    subplot(1,length(files),i); hold on

    % Plot bars
    wake_bar = bar([test_wt(i).wake_avg, test_het(i).wake_avg;0,0]);
    % Set plot settings
    wake_bar(1).FaceColor = color_wt;
    wake_bar(2).FaceColor = color_het;
    wake_bar(1).EdgeColor = 'none';
    wake_bar(2).EdgeColor = 'none';

    % Plot individual readings
    plot(0.85 * ones(length(test_wt(i).wake_all)), test_wt(i).wake_all, Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
    plot(1.15 * ones(length(test_het(i).wake_all)), test_het(i).wake_all, Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

    % Plot standard error
    errorbar([0.85,0.85], test_wt(i).wake_avg, test_wt(i).wake_sem, Color = color_error_bar, lineWidth = line_w_error)
    errorbar([1.15,1.15], test_het(i).wake_avg, test_het(i).wake_sem , Color = color_error_bar, lineWidth = line_w_error)

    xlim([0.5,1.5])
    box off
    ylabel(ylabels{i}, FontName = font_name, FontSize = font_size);
    xticks([0.85,1.15])
    xticklabels({'WT','HET'})
    set(gca, FontName = font_name, FontSize = font_size)
    xtickangle(xlabel_angle)
end

sgtitle('Wake', FontWeight = 'bold', FontName = font_name, FontSize = font_size_title)
hold off

%---------------------------------------------------------------------------------
% NREM
%---------------------------------------------------------------------------------
figure
for i = 1:length(files)
    subplot(1,length(files),i); hold on

    % Plot bars
    nrem_bar = bar([test_wt(i).nrem_avg, test_het(i).nrem_avg;0,0]);
    % Set plot settings
    nrem_bar(1).FaceColor = color_wt;
    nrem_bar(2).FaceColor = color_het;
    nrem_bar(1).EdgeColor = 'none';
    nrem_bar(2).EdgeColor = 'none';

    % Plot individual readings
    plot(0.85 * ones(length(test_wt(i).nrem_all)), test_wt(i).nrem_all, Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
    plot(1.15 * ones(length(test_het(i).nrem_all)), test_het(i).nrem_all, Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

    % Plot standard error
    errorbar([0.85,0.85], test_wt(i).nrem_avg, test_wt(i).nrem_sem, Color = color_error_bar, lineWidth = line_w_error)
    errorbar([1.15,1.15], test_het(i).nrem_avg, test_het(i).nrem_sem , Color = color_error_bar, lineWidth = line_w_error)

    xlim([0.5,1.5])
    box off
    ylabel(ylabels{i}, FontName = font_name, FontSize = font_size);
    xticks([0.85,1.15])
    xticklabels({'WT','HET'})
    set(gca, FontName = font_name, FontSize = font_size)
    xtickangle(xlabel_angle)
end

sgtitle('NREM', FontWeight = 'bold', FontName = font_name, FontSize = font_size_title)
hold off

%---------------------------------------------------------------------------------
% REM
%---------------------------------------------------------------------------------
figure
for i = 1:length(files)
    subplot(1,length(files),i); hold on

    % Plot bars
    rem_bar = bar([test_wt(i).rem_avg, test_het(i).rem_avg;0,0]);
    % Set plot settings
    rem_bar(1).FaceColor = color_wt;
    rem_bar(2).FaceColor = color_het;
    rem_bar(1).EdgeColor = 'none';
    rem_bar(2).EdgeColor = 'none';

    % Plot individual readings
    plot(0.85 * ones(length(test_wt(i).rem_all)), test_wt(i).rem_all, Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)
    plot(1.15 * ones(length(test_het(i).rem_all)), test_het(i).rem_all, Marker = marker,MarkerFaceColor=marker_face_color,MarkerEdgeColor=marker_edge_color,LineStyle = 'none', MarkerSize=marker_size, LineWidth=marker_line_width)

    % Plot standard error
    errorbar([0.85,0.85], test_wt(i).rem_avg, test_wt(i).rem_sem, Color = color_error_bar, lineWidth = line_w_error)
    errorbar([1.15,1.15], test_het(i).rem_avg, test_het(i).rem_sem , Color = color_error_bar, lineWidth = line_w_error)

    xlim([0.5,1.5])
    box off
    ylabel(ylabels{i}, FontName = font_name, FontSize = font_size);
    xticks([0.85,1.15])
    xticklabels({'WT','HET'})
    set(gca, FontName = font_name, FontSize = font_size)
    xtickangle(xlabel_angle)
end

sgtitle('REM', FontWeight = 'bold', FontName = font_name, FontSize = font_size_title)
hold off

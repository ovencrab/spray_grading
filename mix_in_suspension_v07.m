%%% Mix_in_suspension simulation %%%
%%% Author: Jack Evans %%%
%%% Patrick Grant group, Oxford University 2020 %%%
% Simulates spray deposition of a two suspension system where 
% suspension 2 is pumped into suspension 1 and the resulting
% mixture is pumped to the spray head.

clear
clc

%% Input

% Specify coordinates and pattern
load shelves_rotation.mat
load shelves_rotation_defaults.mat

if exist('s1_vol_init','var') == 0
    % Volume of suspension 1 in bottle 1 (ml)
    s1_vol_init = 85;

    % Volume of suspension 2 in bottle 2 (ml)
    s2_vol_b2 = 270;

    % Pump speeds (ml/min)
    p1 = 10;
    p2 = 10;

    % Experiment time
    iteration_total = 1;
    iteration_time = 240*6; % seconds
    time_tot = iteration_total * iteration_time; % seconds

    % Ratio of components
    s1_ratio_binder = 2;
    s1_ratio_cb = 5;
    s1_ratio_am = 100 - s1_ratio_binder - s1_ratio_cb;

    s2_ratio_binder = 2;
    s2_ratio_cb = 0.5;
    s2_ratio_am = 100 - s2_ratio_binder - s2_ratio_cb;

%     save('shelves_rotation_defaults.mat') % delete x_it, x_unique, y_it
%     and y_unique before saving default settings
end

% Overwrite defaults
s2_vol_b2 = 270;
%% Conversions

% Seconds to desired time_step

%conv = 1000; % seconds to millseconds
%time_tot = time_tot*conv;

conv_min_and_s = 60;

p1 = p1 / conv_min_and_s;
p2 = p2 / conv_min_and_s;

mins(:,1) = (1:time_tot) / conv_min_and_s;


%% Pre-allocate

s1_out = zeros(time_tot,1);
s2_out = zeros(time_tot,1);

s1_vol_array = zeros(time_tot,1);
s2_vol_array = zeros(time_tot,1);

out_ratio_am = zeros(time_tot,1);
out_ratio_binder = zeros(time_tot,1);
out_ratio_cb = zeros(time_tot,1);

% Initial volume of suspensions in bottle 1
s1_vol = s1_vol_init;
s2_vol = 0;

% Cumulative volume of suspension 2 pumped into bottle 1
s2_vol_in = 0;

%% Iterate

for i = 1:time_tot
    
    % Proportion of s1 and s2 being pumped from bottle 1 per time_step
    s1_out(i) = p1 * (s1_vol / (s1_vol + s2_vol));
    s2_out(i) = p1 * (s2_vol / (s1_vol + s2_vol));
    
    % Volume of s1_vol left in bottle 1
    if s1_out(i) > s1_vol
        s1_out(i) = s1_vol;
        s1_vol = 0;
    else
        s1_vol = s1_vol - s1_out(i);
    end
    
    % Volume of s2 in bottle 1 considering replenishment from bottle 2
    if s2_vol_in < s2_vol_b2
        if s2_out(i) > s2_vol
            s2_out(i) = s2_vol;
            s2_vol = 0 + p2;
            s2_vol_in = s2_vol_in + p2;
        else
            s2_vol = s2_vol - s2_out(i) + p2;
            s2_vol_in = s2_vol_in + p2;
        end
    else
        if s2_out(i) > s2_vol
            s2_out(i) = s2_vol;
            s2_vol = 0;
        else
            s2_vol = s2_vol - s2_out(i);
        end
    end
    
    % Update volume arrays
    s1_vol_array(i) = s1_vol;
    s2_vol_array(i) = s2_vol;
    
    % Ratio of s1 and s2 out of bottle 1
    s1_out_ratio = s1_out(i) / (s1_out(i) + s2_out(i));
    s2_out_ratio = s2_out(i) / (s1_out(i) + s2_out(i));
    
    % Ratio of components out of bottle 1
    out_ratio_am(i) = s1_ratio_am * s1_out_ratio + s2_ratio_am * s2_out_ratio;
    out_ratio_binder(i) = s1_ratio_binder * s1_out_ratio + s2_ratio_binder * s2_out_ratio;
    out_ratio_cb(i) = s1_ratio_cb * s1_out_ratio + s2_ratio_cb * s2_out_ratio;
end

% Calculate totals and averages

total_vol = s1_vol + s2_vol;
total_out = sum(s1_out) + sum(s2_out);

out_mean_am = mean(out_ratio_am);
out_mean_binder = mean(out_ratio_binder);
out_mean_cb = mean(out_ratio_cb);

%% Time plots

% Formatting
lw = 2;
avg_colour = 'r';
b1_colour = 'm';
b2_colour = 'g';

% Plot 1 -- Volume of suspensions in bottle 1 vs Time
fig_vol = figure('Name', 'Volume in bottle 1', 'Position', [100 100 450 450]);
hold on
plot(mins, s1_vol_array, 'LineWidth', lw, 'DisplayName', 'Suspension 1')
plot(mins, s2_vol_array, 'LineWidth', lw, 'DisplayName', 'Suspension 2')
plot(mins, s1_vol_array+s2_vol_array, '--', 'LineWidth', lw, 'DisplayName', 'Total volume')
legend('Location','best');
set(gca,'FontSize',12)
set(gca,'PlotBoxAspectRatio',[1 1 1])
xlabel('Time (mins)') 
ylabel('Volume (ml)') 
hold off

% Plot 2 -- Composition of sprayed suspension vs Time
fig_out_ratio = figure('Name', 'Sprayed component ratio', 'Position', [100 100 550 450]);

% AM subplot
plot_percent_am = subplot(1,3,1);
plot(mins, out_ratio_am, 'LineWidth', lw, 'DisplayName', 'Transient');
title('AM');
ylim(plot_percent_am,[92.5 98]);
yline(out_mean_am, '--', 'LineWidth', lw, 'Color', avg_colour, 'DisplayName', 'Average');
yline(s1_ratio_am, '-.', 'LineWidth', lw, 'Color', b1_colour, 'DisplayName', 'B1 initial');
yline(s2_ratio_am, '-.', 'LineWidth', lw, 'Color', b2_colour, 'DisplayName', 'B2 initial');
ylabel('Concentration (%)')
set(gca,'FontSize',12)

% Binder subplot
plot_percent_binder = subplot(1,3,2);
plot(mins, out_ratio_binder, 'LineWidth', lw, 'DisplayName', 'Transient');
title('Binder');
ylim(plot_percent_binder,[0 4]);
yline(out_mean_binder, '--', 'LineWidth', lw, 'Color', avg_colour, 'DisplayName', 'Average');
yline(s1_ratio_binder, '-.', 'LineWidth', lw, 'Color', b1_colour, 'DisplayName', 'B1 initial');
yline(s2_ratio_binder, '-.', 'LineWidth', lw, 'Color', b2_colour, 'DisplayName', 'B2 initial');
legend('Location','best');
xlabel('Time (mins)')
set(gca,'FontSize',12)

% CB subplot
plot_percent_cb = subplot(1,3,3);
plot(mins, out_ratio_cb, 'LineWidth', lw, 'DisplayName', 'Transient');
title('CB');
ylim(plot_percent_cb,[0 5.5])
yline(out_mean_cb, '--', 'LineWidth', lw, 'Color', avg_colour, 'DisplayName', 'Average');
yline(s1_ratio_cb, '-.', 'LineWidth', lw, 'Color', b1_colour, 'DisplayName', 'B1 initial');
yline(s2_ratio_cb, '-.', 'LineWidth', lw, 'Color', b2_colour, 'DisplayName', 'B2 initial');
set(gca,'FontSize',12)

%% Spatial plots

% Concatenate iteration coords to total iterations of deposition
x_total = repmat(x_it,1,iteration_total)';
y_total = repmat(y_it,1,iteration_total)';

% Create spline ---------------------------------------------------------

% Number of points desired per iteration
points = length(x_it); 

% Time step for desired number of points
step = 1 * (iteration_time / points) / conv_min_and_s;

% End time of spline
spline_x_end = (iteration_total * iteration_time) / conv_min_and_s;

% Spline distribution
spl_mins = 0+step: step: spline_x_end;

% Spline creation
spl_am = makima(mins, out_ratio_am, spl_mins)';
spl_binder = makima(mins, out_ratio_binder, spl_mins)';
spl_cb = makima(mins, out_ratio_cb, spl_mins)';

spl = [x_total y_total spl_am spl_binder spl_cb];

% Create table of raw data for indexing
table_raw = array2table(spl,'VariableNames',{'x','y','ratio_am','ratio_binder','ratio_cb'});

% Create table for plot data
plot_data = [y_unique' x_unique' zeros(length(x_unique),1)];

for i = 1:length(plot_data)
    rows = (table_raw.y == y_unique(i) & table_raw.x == x_unique(i));
    T = table_raw(rows,:);
    plot_data(i,3) = sum(T.ratio_am)/height(T);
    plot_data(i,4) = sum(T.ratio_binder)/height(T);
    plot_data(i,5) = sum(T.ratio_cb)/height(T);
end

plot_data = array2table(plot_data,'VariableNames',{'y','x','ratio_am','ratio_binder','ratio_cb'});

%%

% Pre-allocate
matrix_am = zeros(max(y_unique),max(x_unique));
matrix_binder = zeros(max(y_unique),max(x_unique));
matrix_cb = zeros(max(y_unique),max(x_unique));

% Build matrix
for i = 1:height(plot_data)
    matrix_am(plot_data.y(i),plot_data.x(i)) = plot_data.ratio_am(i);
    matrix_binder(plot_data.y(i),plot_data.x(i)) = plot_data.ratio_binder(i);
    matrix_cb(plot_data.y(i),plot_data.x(i)) = plot_data.ratio_cb(i);
end

% Round values for display
matrix_am = round(matrix_am, 2);
matrix_binder = round(matrix_binder, 2);
matrix_cb = round(matrix_cb, 2);

% Plot 3 -- Average spatial composition at x vs y coords
fig_heatmaps = figure('Name', 'Component spatial ratio', 'Position', [100 100 800 600]);
t = tiledlayout(2, 2);
nexttile

% AM subplot
h_am = heatmap(matrix_am);
title('AM');
set(gca,'FontSize',12)

% Binder subplot
nexttile
h_binder = heatmap(matrix_binder);
title('Binder');
set(gca,'FontSize',12)

% CB subplot
nexttile
h_cb = heatmap(matrix_cb);
title('CB');
set(gca,'FontSize',12)

t.TileSpacing = 'compact';
t.Padding = 'compact';
xlabel(t,'x dimension');
ylabel(t,'y dimension');

%% Outputs
disp('Suspension 1: ')
disp(['Volume: ',num2str(s1_vol_init),' ml'])
disp(['AM: ',num2str((s1_ratio_am/100)*s1_vol_init/100),' g'])
disp(['Binder: ',num2str((s1_ratio_binder/100)*s1_vol_init/100),' g'])
disp(['CB: ',num2str((s1_ratio_cb/100)*s1_vol_init/100),' g'])
disp(' ')
disp('Suspension 2: ')
disp(['Volume: ',num2str(s2_vol_b2),' ml'])
disp(['AM: ',num2str((s2_ratio_am/100)*s2_vol_b2/100),' g'])
disp(['Binder: ',num2str((s2_ratio_binder/100)*s2_vol_b2/100),' g'])
disp(['CB: ',num2str((s2_ratio_cb/100)*s2_vol_b2/100),' g'])
disp(' ')
disp('Deposition: ')
disp(['Total deposition time: ', num2str(time_tot/conv_min_and_s), ' mins'])
disp(['Final volume in bottle 1: ', num2str(s1_vol + s2_vol), ' ml'])
disp(['Final volume in bottle 2: ', num2str(s2_vol_b2 - s2_vol_in), ' ml'])
disp(['Initial sprayed component ratio: ',...
    num2str(round(out_ratio_am(1),2)),' % AM, ',...
    num2str(round(out_ratio_binder(1),2)),' % Binder, ',...
    num2str(round(out_ratio_cb(1),2)),' % CB'])

disp(['Final sprayed component ratio: ',...
    num2str(round(out_ratio_am(end),2)),' % AM, ',...
    num2str(round(out_ratio_binder(end),2)),' % Binder, ',...
    num2str(round(out_ratio_cb(end),2)),' % CB'])

disp(['Mean component ratio: ',...
    num2str(round(out_mean_am,2)),' % AM, ',...
    num2str(round(out_mean_binder,2)),' % Binder, ',...
    num2str(round(out_mean_cb,2)),' % CB'])

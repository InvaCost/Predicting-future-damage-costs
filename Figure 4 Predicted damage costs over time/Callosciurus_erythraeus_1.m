% --- Raw time data (years of recorded damage costs) ---
time_year = [2001 2002 2003 2005 2006 2007 2008 2009 2010 2011 2012 ...
             2013 2014 2015 2016 2017]';

% Convert years to time since first report (e.g., t = 0 in 2001)
time = time_year - time_year(1);

% Total prediction time horizon (to 2050)
n = 2050 - time_year(1);  % 49 years

% --- Observed cumulative cost data (in USD millions) ---
cost_data = [0.0000227820; 0.0000336500; 0.0000451435; 0.0000679191;
             0.0002142905; 0.0002745281; 0.0003406172; 0.0004383024;
             0.0005151388; 0.0006380282; 0.0007316905; 0.0007748157;
             0.0008334545; 0.0009359811; 0.0009735314; 0.0009913621] * 1e3;

% --- Cost-density function shape parameters ---
s1 = 1;             % shape parameter 1 (threshold-type)
s2 = 0.1;           % shape parameter 2 (steepness of cost escalation)

% --- Precompute constants for cost-density transformation ---
d = (1 - s1) / s2;
B = 1 / (1 + exp(1 / s2 - d));                  % baseline offset
W = (1 + exp(-d)) / (1 - B * (1 + exp(-d)));    % scaling factor

% --- Nonlinear model function definition ---
% Inputs: b1 = [Cmax, gamma, alpha], t = time vector
modelfun1 = @(b1,t) (b1(1) .* W .* ...
    (1 ./ (1 + ((1 - B) / B) .* ...
    exp(-(1 ./ (1 + (b1(2) - 1) .* exp(-b1(3) .* t))) / s2)) - B));

% Initial parameter guesses: [Cmax, gamma, alpha]
beta1 = [1 1 1];

% --- Fit model to cost data using nonlinear least squares ---
mdl1 = fitnlm(time, cost_data, modelfun1, beta1);

% Extract fitted parameters
params = mdl1.Coefficients.Estimate;
Cmax = params(1);    % Maximum cost
gamma = params(2);   % Environmental scaling factor (K / u?)
alp = params(3);     % Intrinsic growth rate

% --- Forecasting over full time horizon ---
t = 0:0.1:n;  % fine time grid for plotting

% Compute fitted cost curve
cost_curve = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B) / B) .* ...
    exp(-(1 ./ (1 + (gamma - 1) .* exp(-alp .* t))) / s2)) - B);

% --- Predict confidence intervals for cost projections ---
tnew = t';
[ynew, ynewci_A] = predict(mdl1, tnew);  % CI from model fit
CI_lower = ynewci_A(:, 1);
CI_upper = ynewci_A(:, 2);

% --- Plot prediction bounds using shaded area ---
hold on;
jbfill(t, cost_curve, CI_upper', [0.3, 0.3, 0.3], 'none', 0, 0.1);

% --- Extract 95% confidence intervals for model parameters ---
ci = round(coefCI(mdl1, 0.05), 3);
Cmax_CI = ci(1, :);
K_CI = ci(2, :);
alp_CI = ci(3, :);

% --- Identify ecological thresholds along the curve ---
term = (1 / s2 - 1) + sqrt((1 / s2 - 1)^2 - 1);
X1 = 1 - s2 * log(term);  % threshold point (rapid rise onset)
t_crit = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X1 - 1));
year_crit = round(time_year(1) + t_crit);
cost_crit = Cmax * W * ...
    (1 / (1 + ((1 - B) / B) * ...
    exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_crit))) / s2)) - B);

% Midpoint (cost = 0.5 * Cmax)
eps2 = 0.5;
X2 = 1 - s2 * log(2 / eps2 - 1);
t_mid = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X2 - 1));
year_mid = round(time_year(1) + t_mid);
cost_mid = Cmax * W * ...
    (1 / (1 + ((1 - B) / B) * ...
    exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_mid))) / s2)) - B);

% Near saturation (cost = 0.9 * Cmax)
eps2 = 0.9;
z_sat = 1 - s2 * log(2 / eps2 - 1);
t_sat = -(1 / alp) * log((1 / (gamma - 1)) * (1 / z_sat - 1));
year_sat = round(time_year(1) + t_sat);
cost_sat = Cmax * W * ...
    (1 / (1 + ((1 - B) / B) * ...
    exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_sat))) / s2)) - B);

% Duration of rapid escalation phase (management window)
management_window = round(t_sat - t_crit, 2);

% --- Cost predictions for last report and 2050 ---
est_cost_at_last_cost_report = round(cost_curve(10 * time(end) + 1), 3);
cost_pred2050 = round(cost_curve(10 * n + 1), 3);
cost_pred2050_upper = round(CI_upper(end), 3);

% Absolute and percentage increases
cost_inc = cost_pred2050 - est_cost_at_last_cost_report;
percent_cost_inc = round((cost_inc / est_cost_at_last_cost_report) * 100, 3);
cost_inc_upper = CI_upper(end) - est_cost_at_last_cost_report;
percent_cost_inc_upper = round((cost_inc_upper / est_cost_at_last_cost_report) * 100, 3);

% --- Plot the results ---
hold on;
plot(t, cost_curve, 'k-', 'linewidth', 2);           % Main fitted curve
plot(time, cost_data, 'kx', 'linewidth', 2);         % Observed data
h4 = plot(time(end), est_cost_at_last_cost_report, 'k^', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);  % Last report
h5 = plot(n, cost_pred2050, 'kp', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);                         % 2050 prediction
h6 = plot(n, CI_upper(end), 'kh', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);                         % Upper bound 2050

% Annotate key ecological points
h1 = plot(t_crit, cost_crit, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(t_mid, cost_mid, 'ks', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 11, 'LineWidth', 2);
h3 = plot(t_sat, cost_sat, 'kd', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Optional legend (commented out)
% legend([h1 h2 h3 h4 h5 h6], {'Threshold point', 'Midpoint', ...
%     'Near saturation point', 'Estimated cost at last cost report', ...
%     'Predicted damage cost in 2050', 'Upper predicted damage cost in 2050'}, ...
%     'Location', 'southeast');
% legend boxoff;

% --- Final plot formatting ---
axis([0 n + 1 0 1.1 * Cmax]);
set(gca, 'fontsize', 18);
title('(a) \it{Callosciurus erythraeus}', 'FontSize', 22, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');
axis square;

% Set x-axis to real years for readability
tick_years = unique([1995:10:2050, 2050]);  % Ensure sorted and unique ticks
set(gca, 'XTick', tick_years - time_year(1), 'XTickLabel', tick_years);

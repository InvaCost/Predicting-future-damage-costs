% --- Time data: years of reported costs ---
time_year = [2000:2017]';  % 18 data points (inclusive)

% Convert to time since first cost report (t = 0 in year 2000)
time = time_year - time_year(1);

% Time horizon up to 2050 (for forecast)
n = 2050 - time_year(1);

% --- Observed cumulative cost data (in USD millions) ---
cost_data = [0.0011227525; 0.0019315128; 0.0027574837; 0.0037918926;
             0.0049553719; 0.0060599913; 0.0072205077; 0.0084654167;
             0.0098312580; 0.0112110614; 0.0124788627; 0.0140081528;
             0.0153328050; 0.0162492165; 0.0168551508; 0.0173421520;
             0.0179523441; 0.0184694342] * 1e3;

% --- Cost-density model shape parameters ---
s1 = 1;         % Threshold form
s2 = 0.1;       % Steepness of cost escalation

% Constants for normalising the cost function
d = (1 - s1) / s2;
B = 1 / (1 + exp(1 / s2 - d));                     % Baseline
W = (1 + exp(-d)) / (1 - B * (1 + exp(-d)));       % Scaling factor

% --- Define the nonlinear model function (cost vs time) ---
modelfun1 = @(b1, t) b1(1) .* W .* ...
    (1 ./ (1 + ((1 - B) / B) .* ...
    exp(-(1 ./ (1 + (b1(2) - 1) .* exp(-b1(3) .* t))) ./ s2)) - B);

% Initial guesses for fitting [Cmax, gamma, alpha]
beta1 = [1 1 1];

% Fit the model to observed data using non-linear regression
mdl1 = fitnlm(time, cost_data, modelfun1, beta1);

% --- Extract fitted parameters ---
params = mdl1.Coefficients.Estimate;
Cmax = params(1);    % Maximum long-term cost
gamma = params(2);   % Environmental scaling factor
alp = params(3);     % Intrinsic growth rate

% --- Forecast the cost curve to 2050 ---
t = 0:0.1:n;
cost_curve = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (gamma - 1).*exp(-alp.*t))) ./ s2)) - B);

% --- Confidence intervals on the forecast ---
tnew = t';
[ynew, ynewci_A] = predict(mdl1, tnew);
CI_lower = ynewci_A(:,1);
CI_upper = ynewci_A(:,2);

% --- Plot the fitted cost curve and confidence region ---
hold on;
jbfill(t, cost_curve, CI_upper', [0.3, 0.3, 0.3], 'none', 0, 0.1);

% Display confidence intervals for model parameters
ci = round(coefCI(mdl1, 0.05), 3);
Cmax_CI = ci(1, :);
K_CI = ci(2, :);
alp_CI = ci(3, :);

% --- Critical ecological/economic thresholds ---
term = (1 / s2 - 1) + sqrt((1 / s2 - 1)^2 - 1);
X1 = 1 - s2 * log(term);  % Threshold density (onset of rapid increase)
t_crit = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X1 - 1));
year_crit = round(time_year(1) + t_crit);
cost_crit = Cmax * W * ...
    (1 / (1 + ((1 - B) / B) * ...
    exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_crit))) / s2)) - B);

% Midpoint (50% of max cost)
eps2 = 0.5;
X2 = 1 - s2 * log(2 / eps2 - 1);
t_mid = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X2 - 1));
year_mid = round(time_year(1) + t_mid);
cost_mid = Cmax * W * ...
    (1 / (1 + ((1 - B) / B) * ...
    exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_mid))) / s2)) - B);

% Near saturation (90% of max cost)
eps2 = 0.9;
z_sat = 1 - s2 * log(2 / eps2 - 1);
t_sat = -(1 / alp) * log((1 / (gamma - 1)) * (1 / z_sat - 1));
year_sat = round(time_year(1) + t_sat);
cost_sat = Cmax * W * ...
    (1 / (1 + ((1 - B) / B) * ...
    exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_sat))) / s2)) - B);

% Time window for management intervention
management_window = round(t_sat - t_crit, 2);

% --- Cost projections ---
est_cost_at_last_cost_report = round(cost_curve(10 * time(end) + 1), 3);
cost_pred2050 = round(cost_curve(10 * n + 1), 3);
cost_pred2050_upper = round(CI_upper(end), 3);

cost_inc = cost_pred2050 - est_cost_at_last_cost_report;
percent_cost_inc = round((cost_inc / est_cost_at_last_cost_report) * 100, 3);
cost_inc_upper = CI_upper(end) - est_cost_at_last_cost_report;
percent_cost_inc_upper = round((cost_inc_upper / est_cost_at_last_cost_report) * 100, 3);

% --- Plot data and key projection points ---
plot(t, cost_curve, 'k-', 'linewidth', 2);              % model curve
plot(time, cost_data, 'kx', 'linewidth', 2);            % raw data
h4 = plot(time(end), est_cost_at_last_cost_report, 'k^', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
h5 = plot(n, cost_pred2050, 'kp', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h6 = plot(n, CI_upper(end), 'kh', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Annotate ecological reference points
h1 = plot(t_crit, cost_crit, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(t_mid, cost_mid, 'ks', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 11, 'LineWidth', 2);
h3 = plot(t_sat, cost_sat, 'kd', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Optional legend
% legend([h1 h2 h3 h4 h5 h6], {'Threshold point', 'Midpoint', ...
%     'Near saturation point', 'Last cost report', ...
%     'Predicted 2050 cost', 'Upper 2050 cost'}, ...
%     'Location', 'southeast');
% legend boxoff;

% --- Formatting ---
axis([0 n+1 0 1.1 * Cmax]);
set(gca, 'fontsize', 18);
title('(c) \it{Myocastor coypus}', 'FontSize', 22, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');
axis square;

% Set real calendar years on x-axis
tick_years = unique([1995:10:2050, 2050]);
set(gca, 'XTick', tick_years - time_year(1), 'XTickLabel', tick_years);

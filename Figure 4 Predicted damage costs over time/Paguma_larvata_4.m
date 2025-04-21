% --- Raw year vector for reported costs ---
time_year = (2000:2017)';  % 18 years of data

% Convert to time since first report (e.g., t = 0 in year 2000)
time = time_year - time_year(1);

% Time horizon until 2050 (for forecast)
n = 2050 - time_year(1);

% --- Cumulative cost data (converted to USD millions) ---
cost_data = [0.0005943984; 0.0015740235; 0.0023347862; 0.0042082156;
             0.0063072763; 0.0083912490; 0.0107959225; 0.0133158593;
             0.0168295962; 0.0207247934; 0.0255526833; 0.0301678621;
             0.0359615432; 0.0406945389; 0.0451999535; 0.0489165414;
             0.0529344219; 0.0566521213] * 1e3;

% --- Cost-density shape parameters ---
s1 = 1;  % type of curve
s2 = 0.1;  % steepness of escalation

% Intermediate constants for model scaling
d = (1 - s1) / s2;
B = 1 / (1 + exp(1 / s2 - d));                      % lower asymptote shift
W = (1 + exp(-d)) / (1 - B * (1 + exp(-d)));        % scaling weight

% --- Define nonlinear model function ---
% Inputs: [Cmax, gamma, alpha], where gamma = K/u0, alpha = growth rate
modelfun1 = @(b1, t) b1(1) .* W .* ...
    (1 ./ (1 + ((1 - B) / B) .* ...
    exp(-(1 ./ (1 + (b1(2) - 1) .* exp(-b1(3) .* t))) ./ s2)) - B);

% Initial guesses for parameters: [Cmax, gamma, alpha]
beta1 = [1 1 1];

% Fit nonlinear regression model
mdl1 = fitnlm(time, cost_data, modelfun1, beta1);

% Extract fitted parameters
params = mdl1.Coefficients.Estimate;
Cmax = params(1);    % long-term maximum cost
gamma = params(2);   % environmental scaling
alp = params(3);     % intrinsic growth rate

% --- Compute fitted cost curve from t = 0 to 2050 ---
t = 0:0.1:n;
cost_curve = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (gamma - 1).*exp(-alp.*t))) ./ s2)) - B);

% --- Predict confidence intervals for each time point ---
tnew = t';
[ynew, ynewci_A] = predict(mdl1, tnew);
CI_lower = ynewci_A(:, 1);
CI_upper = ynewci_A(:, 2);

% --- Plot confidence interval as shaded region ---
hold on;
jbfill(t, cost_curve, CI_upper', [0.3, 0.3, 0.3], 'none', 0, 0.1);  % upper bound only

% --- Extract 95% confidence intervals for model parameters ---
ci = round(coefCI(mdl1, 0.05), 3);
Cmax_CI = ci(1, :);
K_CI = ci(2, :);
alp_CI = ci(3, :);

% --- Compute ecological cost milestones ---
% 1. Threshold time: when costs begin rising rapidly
term = (1 / s2 - 1) + sqrt((1 / s2 - 1)^2 - 1);
X1 = 1 - s2 * log(term);
t_crit = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X1 - 1));
year_crit = round(time_year(1) + t_crit);
cost_crit = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* exp(-(1 ./ (1 + (gamma - 1) .* exp(-alp .* t_crit))) ./ s2)) - B);

% 2. Midpoint: cost = 0.5 * Cmax
eps2 = 0.5;
X2 = 1 - s2 * log(2 / eps2 - 1);
t_mid = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X2 - 1));
year_mid = round(time_year(1) + t_mid);
cost_mid = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* exp(-(1 ./ (1 + (gamma - 1) .* exp(-alp .* t_mid))) ./ s2)) - B);

% 3. Near saturation: cost = 0.9 * Cmax
eps2 = 0.9;
z_sat = 1 - s2 * log(2 / eps2 - 1);
t_sat = -(1 / alp) * log((1 / (gamma - 1)) * (1 / z_sat - 1));
year_sat = round(time_year(1) + t_sat);
cost_sat = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* exp(-(1 ./ (1 + (gamma - 1) .* exp(-alp .* t_sat))) ./ s2)) - B);

% Duration of rapid escalation = management window
management_window = round(t_sat - t_crit, 2);

% --- Forecast costs at present and in 2050 ---
est_cost_at_last_cost_report = round(cost_curve(10 * time(end) + 1), 3);
cost_pred2050 = round(cost_curve(10 * n + 1), 3);
cost_pred2050_upper = round(CI_upper(end), 3);

% Calculate cost increases (absolute and percentage)
cost_inc = cost_pred2050 - est_cost_at_last_cost_report;
percent_cost_inc = round((cost_inc / est_cost_at_last_cost_report) * 100, 3);
cost_inc_upper = CI_upper(end) - est_cost_at_last_cost_report;
percent_cost_inc_upper = round((cost_inc_upper / est_cost_at_last_cost_report) * 100, 3);

% --- Plot cost trajectory and key milestone points ---
plot(t, cost_curve, 'k-', 'linewidth', 2);            % fitted cost curve
plot(time, cost_data, 'kx', 'linewidth', 2);          % observed data
h4 = plot(time(end), est_cost_at_last_cost_report, 'k^', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
h5 = plot(n, cost_pred2050, 'kp', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h6 = plot(n, CI_upper(end), 'kh', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Annotate threshold, midpoint, and saturation points
h1 = plot(t_crit, cost_crit, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(t_mid, cost_mid, 'ks', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 11, 'LineWidth', 2);
h3 = plot(t_sat, cost_sat, 'kd', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% --- Final formatting ---
axis([0 n+1 0 1.1 * Cmax]);
set(gca, 'fontsize', 18);
title('(d) \it{Paguma larvata}', 'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');
axis square;

% Real year x-axis labels
tick_years = unique([1995:10:2050, 2050]);
set(gca, 'XTick', tick_years - time_year(1), 'XTickLabel', tick_years);

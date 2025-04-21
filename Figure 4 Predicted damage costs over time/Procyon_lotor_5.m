% --- Raw time data (years of reported costs) ---
time_year = (2000:2017)';  % 18 years

% Convert to time since first report (t = 0 in 2000)
time = time_year - time_year(1);

% Forecast horizon until 2050
n = 2050 - time_year(1);  % 50 years

% --- Cumulative cost data (USD millions) ---
cost_data = [0.0004755187; 0.0008855944; 0.0017333013; 0.0026412825;
             0.0041885900; 0.0059537035; 0.0076683403; 0.0097866935;
             0.0119456040; 0.0153523752; 0.0198601133; 0.0250897390;
             0.0295453875; 0.0332002521; 0.0364644787; 0.0394035735;
             0.0425577974; 0.0454731157] * 1e3;

% --- Cost-density function shape parameters ---
s1 = 1;
s2 = 0.1;

% Intermediate constants for cost scaling
d = (1 - s1) / s2;
B = 1 / (1 + exp(1 / s2 - d));
W = (1 + exp(-d)) / (1 - B * (1 + exp(-d)));

% --- Define nonlinear model function ---
modelfun1 = @(b1, t) b1(1) .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (b1(2) - 1) .* exp(-b1(3) .* t))) ./ s2)) - B);

% Initial guesses for [Cmax, gamma, alpha]
beta1 = [1 1 1];
mdl1 = fitnlm(time, cost_data, modelfun1, beta1);  % fit nonlinear model

% Extract estimated parameters
params = mdl1.Coefficients.Estimate;
Cmax = params(1);    % long-term cost saturation
gamma = params(2);   % environmental scaling factor
alp = params(3);     % intrinsic growth rate

% --- Fitted cost curve over full time horizon ---
t = 0:0.1:n;
cost_curve = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (gamma - 1).*exp(-alp.*t))) ./ s2)) - B);

% --- Confidence intervals for predicted cost curve ---
tnew = t';
[ynew, ynewci_A] = predict(mdl1, tnew);
CI_lower = ynewci_A(:, 1);
CI_upper = ynewci_A(:, 2);

% --- Plot shaded confidence region ---
hold on;
jbfill(t, cost_curve, CI_upper', [0.3, 0.3, 0.3], 'none', 0, 0.1);

% --- Extract parameter confidence intervals ---
ci = round(coefCI(mdl1, 0.05), 3);
Cmax_CI = ci(1, :);
K_CI = ci(2, :);
alp_CI = ci(3, :);

% --- Ecological milestones ---
% Threshold time: where cost starts increasing sharply
term = (1 / s2 - 1) + sqrt((1 / s2 - 1)^2 - 1);
X1 = 1 - s2 * log(term);
t_crit = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X1 - 1));
year_crit = round(time_year(1) + t_crit);
cost_crit = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (gamma - 1) .* exp(-alp .* t_crit))) ./ s2)) - B);

% Midpoint time: where cost = 0.5 * Cmax
eps2 = 0.5;
X2 = 1 - s2 * log(2 / eps2 - 1);
t_mid = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X2 - 1));
year_mid = round(time_year(1) + t_mid);
cost_mid = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (gamma - 1) .* exp(-alp .* t_mid))) ./ s2)) - B);

% Saturation time: where cost = 0.9 * Cmax
eps2 = 0.9;
z_sat = 1 - s2 * log(2 / eps2 - 1);
t_sat = -(1 / alp) * log((1 / (gamma - 1)) * (1 / z_sat - 1));
year_sat = round(time_year(1) + t_sat);
cost_sat = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (gamma - 1) .* exp(-alp .* t_sat))) ./ s2)) - B);

% Management window = period between threshold and saturation
management_window = round(t_sat - t_crit, 2);

% --- Projected cost metrics ---
est_cost_at_last_cost_report = round(cost_curve(10 * time(end) + 1), 3);
cost_pred2050 = round(cost_curve(10 * n + 1), 3);
cost_pred2050_upper = round(CI_upper(end), 3);

% Absolute and percentage increase
cost_inc = cost_pred2050 - est_cost_at_last_cost_report;
percent_cost_inc = round((cost_inc / est_cost_at_last_cost_report) * 100, 3);
cost_inc_upper = CI_upper(end) - est_cost_at_last_cost_report;
percent_cost_inc_upper = round((cost_inc_upper / est_cost_at_last_cost_report) * 100, 3);

% --- Plot fitted curve and key points ---
plot(t, cost_curve, 'k-', 'linewidth', 2);               % fitted trajectory
plot(time, cost_data, 'kx', 'linewidth', 2);             % actual data
h4 = plot(time(end), est_cost_at_last_cost_report, 'k^', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
h5 = plot(n, cost_pred2050, 'kp', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h6 = plot(n, CI_upper(end), 'kh', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Critical ecological points
h1 = plot(t_crit, cost_crit, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(t_mid, cost_mid, 'ks', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 11, 'LineWidth', 2);
h3 = plot(t_sat, cost_sat, 'kd', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% --- Final axis and label formatting ---
axis([0 n + 1 0 1.1 * Cmax]);
set(gca, 'fontsize', 18);
title('(e) \it{Procyon lotor}', 'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');
axis square;

% Set x-axis to calendar years
tick_years = unique([1995:10:2050, 2050]);
set(gca, 'XTick', tick_years - time_year(1), 'XTickLabel', tick_years);

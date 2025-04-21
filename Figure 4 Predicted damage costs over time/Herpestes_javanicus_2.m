% --- Time of reported costs (in years) ---
time_year = [2000 2001 2002 2003 2008 2010 2011 2016 2017]';

% Convert to time since first report (t = 0 in 2000)
time = time_year - time_year(1);

% Time span to 2050 (for prediction)
n = 2050 - time_year(1);  % 50 years

% --- Cost data (cumulative), converted to USD millions ---
cost_data = [0.0001056708; 0.0001398438; 0.0001507118; 0.0001622053;
             0.0001732201; 0.0001988323; 0.0002261410; 0.0002355286;
             0.0002444439] * 1e3;

% --- Cost-density function parameters ---
s1 = 1;         % threshold-type cost curve
s2 = 0.1;       % steepness of cost escalation

% Compute constants for logistic cost-density function
d = (1 - s1) / s2;
B = 1 / (1 + exp(1 / s2 - d));                      % offset
W = (1 + exp(-d)) / (1 - B * (1 + exp(-d)));        % scaling

% --- Define nonlinear model function: [Cmax, gamma, alpha] ---
modelfun1 = @(b1,t) (b1(1) .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* ...
    exp(-(1 ./ (1 + (b1(2) - 1) .* exp(-b1(3) .* t))) / s2)) - B));

beta1 = [1 1 1];  % initial parameter estimates
mdl1 = fitnlm(time, cost_data, modelfun1, beta1);  % fit model

% --- Extract estimated parameters ---
params = mdl1.Coefficients.Estimate;
Cmax = params(1);    % maximum cost
gamma = params(2);   % env. scaling factor
alp = params(3);     % growth rate

% --- Generate fitted curve over full time range to 2050 ---
t = 0:0.1:n;
cost_curve = Cmax .* W .* ...
    (1 ./ (1 + ((1 - B)/B) .* exp(-(1 ./ (1 + (gamma - 1).*exp(-alp.*t))) / s2)) - B);

% --- Confidence intervals around fitted predictions ---
tnew = t';
[ynew, ynewci_A] = predict(mdl1, tnew);
CI_lower = ynewci_A(:,1);
CI_upper = ynewci_A(:,2);

% --- Plot shaded confidence region ---
hold on;
jbfill(t, cost_curve, CI_upper', [0.3, 0.3, 0.3], 'none', 0, 0.1);

% --- Parameter confidence intervals ---
ci = round(coefCI(mdl1, 0.05), 3);
Cmax_CI = ci(1, :);
K_CI = ci(2, :);
alp_CI = ci(3, :);

% --- Ecological thresholds: rapid rise, midpoint, saturation ---
term = (1/s2 - 1) + sqrt((1/s2 - 1)^2 - 1);
X1 = 1 - s2 * log(term);  % threshold density
t_crit = -(1/alp) * log((1/(gamma - 1)) * (1/X1 - 1));
year_crit = round(time_year(1) + t_crit);
cost_crit = Cmax * W * ...
    (1 / (1 + ((1 - B)/B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_crit))) / s2)) - B);

% Midpoint (50% of Cmax)
eps2 = 0.5;
X2 = 1 - s2 * log(2/eps2 - 1);
t_mid = -(1/alp) * log((1/(gamma - 1)) * (1/X2 - 1));
year_mid = round(time_year(1) + t_mid);
cost_mid = Cmax * W * ...
    (1 / (1 + ((1 - B)/B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_mid))) / s2)) - B);

% Saturation point (90% of Cmax)
eps2 = 0.9;
z_sat = 1 - s2 * log(2/eps2 - 1);
t_sat = -(1/alp) * log((1/(gamma - 1)) * (1/z_sat - 1));
year_sat = round(time_year(1) + t_sat);
cost_sat = Cmax * W * ...
    (1 / (1 + ((1 - B)/B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_sat))) / s2)) - B);

% Duration of management window (saturation - threshold)
management_window = round(t_sat - t_crit, 2);

% --- Forecasted costs ---
est_cost_at_last_cost_report = round(cost_curve(10 * time(end) + 1), 3);
cost_pred2050 = round(cost_curve(10 * n + 1), 3);
cost_pred2050_upper = round(CI_upper(end), 3);

cost_inc = cost_pred2050 - est_cost_at_last_cost_report;
percent_cost_inc = round((cost_inc / est_cost_at_last_cost_report) * 100, 3);
cost_inc_upper = CI_upper(end) - est_cost_at_last_cost_report;
percent_cost_inc_upper = round((cost_inc_upper / est_cost_at_last_cost_report) * 100, 3);

% --- Plotting the curve and key points ---
hold on;
plot(t, cost_curve, 'k-', 'linewidth', 2);         % fitted curve
plot(time, cost_data, 'kx', 'linewidth', 2);       % actual data
h4 = plot(time(end), est_cost_at_last_cost_report, 'k^', ...
    'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
h5 = plot(n, cost_pred2050, 'kp', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h6 = plot(n, CI_upper(end), 'kh', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Annotate critical points on curve
h1 = plot(t_crit, cost_crit, 'ko', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(t_mid, cost_mid, 'ks', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 11, 'LineWidth', 2);
h3 = plot(t_sat, cost_sat, 'kd', 'MarkerFaceColor', 'k', ...
    'MarkerSize', 10, 'LineWidth', 2);

% Optional legend (commented out)
% legend([h1 h2 h3 h4 h5 h6], {'Threshold point', 'Midpoint', ...
%     'Near saturation point', 'Last reported cost', ...
%     'Predicted 2050 cost', 'Upper 2050 cost'}, ...
%     'Location', 'southeast');
% legend boxoff;

% --- Formatting and axis settings ---
axis([0 n+1 0 1.5 * Cmax]);
set(gca, 'fontsize', 18);
title('(b) \it{Herpestes javanicus}', 'FontSize', 22, ...
    'FontWeight', 'bold', 'Interpreter', 'latex');
axis square;

% Label x-axis in actual years
tick_years = unique([1995:10:2050, 2050]);
set(gca, 'XTick', tick_years - time_year(1), 'XTickLabel', tick_years);

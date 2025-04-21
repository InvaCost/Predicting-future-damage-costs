% Time vector for observed cost data (t = 0 to 9 years)
time = [0; 1; 2; 3; 4; 5; 6; 7; 8; 9];

% Number of years to simulate prediction
n = 100;

% Randomly generated cumulative cost data (heuristic example)
cost_data = [0.706; 0.7378; 1.0147; 1.0609; 1.158; 1.9815; 2.6763; 2.9934; 3.9436; 3.978];

% Cost-density function shape parameters (s1 = 1, high-threshold; s2 = 0.1, sharp escalation)
s1 = 1;
s2 = 0.1;

% Precompute terms for cost-density function
d = (1 - s1) / s2;
B = 1 / (1 + exp(1 / s2 - d));  % Inflection point of cost-density curve
W = (1 + exp(-d)) / (1 - B * (1 + exp(-d)));  % Normalisation term for scaling

% Nonlinear model function combining logistic growth and cost-density relationship
modelfun1 = @(b1, t) (b1(1) .* W .* (1 ./ (1 + ((1 - B) ./ B) .* exp(-(1 ./ (1 + (b1(2) - 1) .* exp(-b1(3) .* t))) ./ s2)) - B));
beta1 = [1 1 1];  % Initial parameter estimates [Cmax, gamma, alpha]

% Fit nonlinear regression model to cost data
mdl1 = fitnlm(time, cost_data, modelfun1, beta1);

% Extract estimated parameters from fitted model
params = mdl1.Coefficients.Estimate;
Cmax = params(1);   % Long-term cost ceiling
gamma = params(2);  % Environmental scaling factor = K / u0
alp = params(3);    % Intrinsic population growth rate

% Predict costs over a finer time grid (0 to 100 years)
t = 0:0.1:n;
cost_curve = Cmax * W * (1 ./ (1 + ((1 - B) / B) .* exp(-(1 ./ (1 + (gamma - 1) * exp(-alp * t))) ./ s2)) - B);

% --- Identify key thresholds along cost curve --- %
% 1. Threshold point: point of rapid cost increase
term = (1 / s2 - 1) + sqrt((1 / s2 - 1)^2 - 1);
X1 = 1 - s2 * log(term);
t_crit1 = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X1 - 1));
cost_crit1 = Cmax * W * (1 / (1 + ((1 - B) / B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_crit1))) / s2)) - B);

% 2. Midpoint: when cost reaches 50% of Cmax
eps2 = 0.5;
X2 = 1 - s2 * log(2 / eps2 - 1);
t_mid1 = -(1 / alp) * log((1 / (gamma - 1)) * (1 / X2 - 1));
cost_mid1 = Cmax * W * (1 / (1 + ((1 - B) / B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_mid1))) / s2)) - B);

% 3. Near-saturation point: cost reaches 90% of Cmax
eps2 = 0.9;
z_sat = 1 - s2 * log(2 / eps2 - 1);
t_sat1 = -(1 / alp) * log((1 / (gamma - 1)) * (1 / z_sat - 1));
cost_sat1 = Cmax * W * (1 / (1 + ((1 - B) / B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp * t_sat1))) / s2)) - B);

% Duration of rapid cost escalation (management window)
management_window = round(t_sat1 - t_crit1, 2);

% Predicted final cost and future projection
est_cost_at_last_cost_report = round(cost_curve(10 * time(end) + 1), 3);
cost_pred2050 = round(cost_curve(10 * n + 1), 3);
cost_inc = cost_pred2050 - est_cost_at_last_cost_report;
percent_cost_inc = round((cost_inc / est_cost_at_last_cost_report) * 100, 3);

% --- Plotting --- %
hold on;
plot(t, cost_curve, 'k-', 'linewidth', 2);                  % Solid curve for original model
plot(time, cost_data, 'kx', 'linewidth', 2);                % Observed data
plot(n, cost_curve(10 * n + 1), 'k^', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2); % Predicted 2050

% Annotate threshold, midpoint, and saturation points
h1 = plot(t_crit1, cost_crit1, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
h2 = plot(t_mid1, cost_mid1, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 11, 'LineWidth', 2);
h3 = plot(t_sat1, cost_sat1, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
legend([h1 h2 h3], {'Threshold point', 'Midpoint', 'Near saturation point'}, 'Location', 'southeast');
legend boxoff;

axis([0 60 0 12]);
set(gca, 'fontsize', 18);
axis square;

% --- Second scenario with reduced alpha (slower growth) --- %
alp2 = 0.5 * alp;

% Repeat same threshold calculations for reduced growth rate
t_crit2 = -(1 / alp2) * log((1 / (gamma - 1)) * (1 / X1 - 1));
cost_crit2 = Cmax * W * (1 / (1 + ((1 - B) / B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp2 * t_crit2))) / s2)) - B);

t_mid2 = -(1 / alp2) * log((1 / (gamma - 1)) * (1 / X2 - 1));
cost_mid2 = Cmax * W * (1 / (1 + ((1 - B) / B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp2 * t_mid2))) / s2)) - B);

t_sat2 = -(1 / alp2) * log((1 / (gamma - 1)) * (1 / z_sat - 1));
cost_sat2 = Cmax * W * (1 / (1 + ((1 - B) / B) * exp(-(1 / (1 + (gamma - 1) * exp(-alp2 * t_sat2))) / s2)) - B);

cost_curve2 = Cmax * W * (1 ./ (1 + ((1 - B) / B) .* exp(-(1 ./ (1 + (gamma - 1) * exp(-alp2 * t))) ./ s2)) - B);

% Plot second scenario (slower growth) as dashed curve
plot(t, cost_curve2, 'k--', 'linewidth', 2);
plot(t_crit2, cost_crit2, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);
plot(t_mid2, cost_mid2, 'ks', 'MarkerFaceColor', 'k', 'MarkerSize', 11, 'LineWidth', 2);
plot(t_sat2, cost_sat2, 'kd', 'MarkerFaceColor', 'k', 'MarkerSize', 10, 'LineWidth', 2);

management_window2 = round(t_sat2 - t_crit2, 2);

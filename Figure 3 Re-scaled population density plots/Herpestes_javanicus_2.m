% Define the species being modelled
species = {'Herpestes javanicus'};

% Logistic growth model parameters (? and ?) and their 95% confidence intervals
gamma_vals = [1.15];                   % Environmental scaling factor ? = K/u?
gamma_conf = [1.07, 1.24];             % Confidence interval for ?

alpha_vals = [0.085];                  % Intrinsic growth rate ? (year?¹)
alpha_conf = [-0.032, 0.201];          % Confidence interval for ?

% Define the time range from 30 years before to 20 years after first cost report
t = -30:0.01:20;

% Begin plotting
figure;
for i = 1:length(species)
    gamma = gamma_vals(i);
    alpha = alpha_vals(i);
    
    % --- Central logistic curve ---
    % Equation (5): z(t) = 1 / [1 + (? ? 1) * exp(??t)]
    z = 1 ./ (1 + (gamma - 1) .* exp(-alpha .* t));
    
    % --- Confidence bounds ---
    gamma_low = gamma_conf(i, 1);
    gamma_high = gamma_conf(i, 2);
    alpha_low = alpha_conf(i, 1);
    alpha_high = alpha_conf(i, 2);
    
    % Set a floor to ensure non-negative ? (for meaningful biological interpretation)
    alpha_low = max(alpha_low, 0);

    % Lower and upper bounds for z(t) across time
    z_low = 1 ./ (1 + (gamma_low - 1) .* exp(-alpha_low .* t));
    z_high = 1 ./ (1 + (gamma_high - 1) .* exp(-alpha_high .* t));

    % --- Inflection point ---
    % Where z(t) = 0.5 ? t_inflection = (1/?) * log(? ? 1)
    tinf = (1 / alpha) * log(gamma - 1);  % For H. javanicus, this is < 0

    % Begin plotting
    hold on;

    % --- Confidence region (shaded band) ---
    c = [0.3 0.3 0.3];  % grey
    fill([t, fliplr(t)], [z_low, fliplr(z_high)], c, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.1);  % semi-transparent fill
    
    % --- Plot central logistic curve ---
    plot(t, z, 'k-', 'linewidth', 2.5);  % central re-scaled density curve
    
    % --- Annotate key points ---
    % (1) Initial re-scaled density at t = 0 (i.e., z(0) = 1/?)
    h1 = plot(0, 1 / gamma, 'r+', 'MarkerFaceColor', 'r', 'MarkerSize', 13, 'LineWidth', 2.5);

    % (2) Inflection point (where growth rate is maximal)
    h2 = plot(tinf, 0.5, 'k*', 'MarkerFaceColor', 'k', 'MarkerSize', 13, 'LineWidth', 2.5);

    % --- Add legend with LaTeX formatting ---
    legend([h1 h2], {'Initial re-scaled density $z(0)=0.868$', ...
                     'Inflection point $z(-22.32)=0.5$'}, ...
           'Location', 'northwest', 'Interpreter', 'latex');
    legend boxoff;

    % --- Axis and formatting ---
    axis([t(1) t(end) 0 1]);              % y-axis from 0 to 1
    set(gca, 'YTick', [], 'XTick', []);  % Suppress tick marks for clean look
    set(gca, 'fontsize', 18);            % Increase font size
    title(['(b) \it{' species{i} '}'], ...
           'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');
    axis square;                         % Square aspect ratio
    
    hold off;
end

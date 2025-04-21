% Define the species being plotted
species = {'Callosciurus erythraeus'};

% Model parameters (? and ?) and their 95% confidence intervals
gamma_vals = [1.98];                   % Environmental scaling factor ? = K/u?
gamma_conf = [1.80, 2.16];             % Confidence interval for ?

alpha_vals = [0.217];                  % Intrinsic population growth rate ?
alpha_conf = [0.185, 0.249];           % Confidence interval for ?

% Define a wide time range centred at t = 0 (first reported cost), spanning 50 years total
t = -30:0.01:20;

figure;
for i = 1:length(species)
    gamma = gamma_vals(i);
    alpha = alpha_vals(i);
    
    % --- Central population growth curve z(t) ---
    % Equation (5): z(t) = 1 / [1 + (? ? 1) * exp(??t)]
    z = 1 ./ (1 + (gamma - 1) .* exp(-alpha .* t));
    
    % --- Confidence region: bounds on z(t) based on parameter uncertainty ---
    gamma_low = gamma_conf(i, 1);
    gamma_high = gamma_conf(i, 2);
    alpha_low = alpha_conf(i, 1);
    alpha_high = alpha_conf(i, 2);

    % Ensure alpha is not negative (for robustness)
    alpha_low = max(alpha_low, 0);

    % Lower and upper bound trajectories of z(t)
    z_low = 1 ./ (1 + (gamma_low - 1) .* exp(-alpha_low .* t));
    z_high = 1 ./ (1 + (gamma_high - 1) .* exp(-alpha_high .* t));

    % --- Inflection point ---
    % Occurs where z = 0.5 ? t_inflection = (1/?) * log(? - 1)
    tinf = (1 / alpha) * log(gamma - 1);

    % Begin plotting
    hold on;
    
    % --- Confidence region (shaded area between z_low and z_high) ---
    c = [0.3 0.3 0.3];  % Grey fill colour
    fill([t, fliplr(t)], [z_low, fliplr(z_high)], c, ...
         'EdgeColor', 'none', 'FaceAlpha', 0.1);  % Transparent shaded band
    
    % --- Central curve ---
    plot(t, z, 'k-', 'linewidth', 2.5);  % Main logistic growth curve
    
    % --- Key points ---
    % (1) Initial re-scaled density at t = 0 is z(0) = 1/?
    h1 = plot(0, 1 / gamma, 'r+', 'MarkerFaceColor', 'r', 'MarkerSize', 13, 'LineWidth', 2.5);
    
    % (2) Inflection point where z = 0.5
    h2 = plot(tinf, 0.5, 'k*', 'MarkerFaceColor', 'k', 'MarkerSize', 13, 'LineWidth', 2.5);

    % --- Legend using LaTeX interpreter ---
    legend([h1 h2], {'Initial re-scaled density $z(0)=0.505$', ...
                     'Inflection point $z(-0.09)=0.5$'}, ...
           'Location', 'northwest', 'Interpreter', 'latex');
    legend boxoff;

    % --- Axis formatting ---
    axis([t(1) t(end) 0 1]);            % Set limits
    set(gca, 'YTick', 0:0.2:1);         % y-axis ticks
    set(gca, 'XTick', []);             % x-axis ticks removed for clarity
    set(gca, 'fontsize', 18);          % axis font size
    title(['(a) \it{' species{i} '}'], 'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');
    axis square;                       % Make plot square

    hold off;
end

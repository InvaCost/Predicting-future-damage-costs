% Define the species to be plotted
species = {'Myocastor coypus'};

% Parameters for the logistic growth model and their confidence intervals
gamma_vals = [1.52];                   % Environmental scaling factor ? = K/u?
gamma_conf = [1.50, 1.55];             % Confidence interval for ?

alpha_vals = [0.168];                  % Intrinsic growth rate ? (year?¹)
alpha_conf = [0.156, 0.179];           % Confidence interval for ?

% Define the time range: 30 years before to 20 years after t = 0 (first cost report)
t = -30:0.01:20;

% Begin figure
figure;
for i = 1:length(species)
    gamma = gamma_vals(i);
    alpha = alpha_vals(i);

    % --- Central logistic growth curve z(t) ---
    % Equation: z(t) = 1 / [1 + (? ? 1) * exp(??t)]
    z = 1 ./ (1 + (gamma - 1) .* exp(-alpha .* t));
    
    % --- Confidence bounds ---
    gamma_low = gamma_conf(i, 1);
    gamma_high = gamma_conf(i, 2);
    alpha_low = alpha_conf(i, 1);
    alpha_high = alpha_conf(i, 2);
    
    % Prevent negative growth rate (for robustness)
    alpha_low = max(alpha_low, 0);

    % Compute upper and lower confidence region bounds for z(t)
    z_low = 1 ./ (1 + (gamma_low - 1) .* exp(-alpha_low .* t));
    z_high = 1 ./ (1 + (gamma_high - 1) .* exp(-alpha_high .* t));

    % --- Inflection point (where z(t) = 0.5) ---
    % tinf = (1/?) * log(? - 1)
    tinf = (1 / alpha) * log(gamma - 1);

    % --- Begin plot configuration ---
    hold on;

    % --- Confidence region (shaded area) ---
    c = [0.3 0.3 0.3];  % Grey fill
    fill([t, fliplr(t)], [z_low, fliplr(z_high)], c, ...
         'EdgeColor', 'none', 'FaceAlpha', 0.1);  % Transparent bounds
    
    % --- Central logistic curve ---
    plot(t, z, 'k-', 'linewidth', 2.5);  % Solid black curve for mean trajectory
    
    % --- Annotate key points ---
    % (1) Initial re-scaled density at t = 0: z(0) = 1/? = 0.656
    h1 = plot(0, 1 / gamma, 'r+', 'MarkerFaceColor', 'r', ...
              'MarkerSize', 13, 'LineWidth', 2.5);

    % (2) Inflection point at z = 0.5
    h2 = plot(tinf, 0.5, 'k*', 'MarkerFaceColor', 'k', ...
              'MarkerSize', 13, 'LineWidth', 2.5);

    % --- Legend with LaTeX formatting ---
    legend([h1 h2], {'Initial re-scaled density $z(0)=0.656$', ...
                     'Inflection point $z(-3.89)=0.5$'}, ...
           'Location', 'northwest', 'Interpreter', 'latex');
    legend boxoff;

    % --- Axis and styling ---
    axis([t(1) t(end) 0 1]);                     % Fix plot limits
    set(gca, 'YTick', [], 'XTick', t(1):10:t(end)); % Custom x-tick spacing
    set(gca, 'fontsize', 18);                    % Increase font size
    title(['(c) \it{' species{i} '}'], ...
           'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');
    axis square;                                % Maintain square aspect ratio

    hold off;
end

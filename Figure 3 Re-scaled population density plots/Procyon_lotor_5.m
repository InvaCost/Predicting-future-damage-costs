% Define the species being plotted
species = {'Procyon lotor'};

% Logistic growth model parameters and their 95% confidence intervals
gamma_vals = [2.09];                   % Environmental scaling factor ? = K/u?
gamma_conf = [1.98, 2.19];             % Confidence interval for ?

alpha_vals = [0.172];                  % Intrinsic population growth rate ? (year?¹)
alpha_conf = [0.156, 0.188];           % Confidence interval for ?

% Define a broad time range around the first reported cost (t = 0)
t = -30:0.01:20;

% Start the figure
figure;
for i = 1:length(species)
    gamma = gamma_vals(i);
    alpha = alpha_vals(i);

    % --- Central logistic population growth curve ---
    % Equation: z(t) = 1 / [1 + (? ? 1) * exp(??t)]
    z = 1 ./ (1 + (gamma - 1) .* exp(-alpha .* t));

    % --- Confidence intervals ---
    gamma_low = gamma_conf(i, 1);
    gamma_high = gamma_conf(i, 2);
    alpha_low = alpha_conf(i, 1);
    alpha_high = alpha_conf(i, 2);

    % Ensure ? stays non-negative for biological realism
    alpha_low = max(alpha_low, 0);

    % Calculate lower and upper bounds for z(t)
    z_low = 1 ./ (1 + (gamma_low - 1) .* exp(-alpha_low .* t));
    z_high = 1 ./ (1 + (gamma_high - 1) .* exp(-alpha_high .* t));

    % --- Inflection point ---
    % Occurs when z(t) = 0.5 ? t_inf = (1/?) * log(? ? 1)
    tinf = (1 / alpha) * log(gamma - 1);

    % --- Plotting starts ---
    hold on;

    % --- Shaded confidence region between z_low and z_high ---
    c = [0.3 0.3 0.3];  % grey fill colour
    fill([t, fliplr(t)], [z_low, fliplr(z_high)], c, ...
         'EdgeColor', 'none', 'FaceAlpha', 0.1);  % transparent shading

    % --- Central growth curve ---
    plot(t, z, 'k-', 'linewidth', 2.5);  % main z(t) trajectory

    % --- Annotate key points ---
    % (1) Initial re-scaled density at t = 0: z(0) = 1/? = 0.479
    h1 = plot(0, 1 / gamma, 'r+', 'MarkerFaceColor', 'r', ...
              'MarkerSize', 13, 'LineWidth', 2.5);

    % (2) Inflection point at z = 0.5, t ? 0.50
    h2 = plot(tinf, 0.5, 'k*', 'MarkerFaceColor', 'k', ...
              'MarkerSize', 13, 'LineWidth', 2.5);

    % --- Legend using LaTeX formatting ---
    legend([h1 h2], {'Initial re-scaled density $z(0)=0.479$', ...
                     'Inflection point $z(0.50)=0.5$'}, ...
           'Location', 'northwest', 'Interpreter', 'latex');
    legend boxoff;

    % --- Axis formatting and labels ---
    axis([t(1) t(end) 0 1]);                  % Set axis limits
    set(gca, 'YTick', [], 'XTick', t(1):10:t(end));  % Show only x-ticks
    set(gca, 'fontsize', 18);                % Increase font size
    title(['(e) \it{' species{i} '}'], ...
           'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');
    axis square;                             % Keep plot square

    hold off;
end

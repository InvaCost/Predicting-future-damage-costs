% Define the species to be plotted
species = {'Paguma larvata'};

% Logistic model parameters and their 95% confidence intervals
gamma_vals = [1.93];                    % Environmental scaling factor ? = K/u?
gamma_conf = [1.88, 1.97];              % Confidence interval for ?

alpha_vals = [0.151];                   % Intrinsic growth rate ? (year?¹)
alpha_conf = [0.142, 0.161];            % Confidence interval for ?

% Define a wide time range: 30 years before to 20 years after first cost report
t = -30:0.01:20;

% Start figure
figure;
for i = 1:length(species)
    gamma = gamma_vals(i);
    alpha = alpha_vals(i);

    % --- Central logistic growth curve z(t) ---
    % Equation: z(t) = 1 / [1 + (? ? 1) * exp(??t)]
    z = 1 ./ (1 + (gamma - 1) .* exp(-alpha .* t));

    % --- Confidence bounds based on ? and ? uncertainty ---
    gamma_low = gamma_conf(i, 1);
    gamma_high = gamma_conf(i, 2);
    alpha_low = alpha_conf(i, 1);
    alpha_high = alpha_conf(i, 2);

    % Ensure non-negative lower bound for ?
    alpha_low = max(alpha_low, 0);

    % Lower and upper bound curves for z(t)
    z_low = 1 ./ (1 + (gamma_low - 1) .* exp(-alpha_low .* t));
    z_high = 1 ./ (1 + (gamma_high - 1) .* exp(-alpha_high .* t));

    % --- Inflection point ---
    % Occurs where z(t) = 0.5 ? t_inf = (1/?) * log(? ? 1)
    tinf = (1 / alpha) * log(gamma - 1);

    % --- Plotting configuration ---
    hold on;

    % --- Confidence region (shaded area between z_low and z_high) ---
    c = [0.3 0.3 0.3];  % grey colour
    fill([t, fliplr(t)], [z_low, fliplr(z_high)], c, ...
         'EdgeColor', 'none', 'FaceAlpha', 0.1);  % light transparent shading

    % --- Plot the central logistic curve ---
    plot(t, z, 'k-', 'linewidth', 2.5);  % central population growth trajectory

    % --- Annotate key reference points ---
    % (1) Initial re-scaled population density at t = 0: z(0) = 1/? = 0.519
    h1 = plot(0, 1 / gamma, 'r+', 'MarkerFaceColor', 'r', ...
              'MarkerSize', 13, 'LineWidth', 2.5);

    % (2) Inflection point at z = 0.5
    h2 = plot(tinf, 0.5, 'k*', 'MarkerFaceColor', 'k', ...
              'MarkerSize', 13, 'LineWidth', 2.5);

    % --- Legend with LaTeX formatting for z-values ---
    legend([h1 h2], {'Initial re-scaled density $z(0)=0.519$', ...
                     'Inflection point $z(-0.48)=0.5$'}, ...
           'Location', 'northwest', 'Interpreter', 'latex');
    legend boxoff;

    % --- Axis configuration ---
    axis([t(1) t(end) 0 1]);                  % fix plot limits
    set(gca, 'YTick', 0:0.2:1, 'XTick', t(1):10:t(end));  % tick marks
    set(gca, 'fontsize', 18);                % axis label font size
    title(['(d) \it{' species{i} '}'], ...
           'FontSize', 22, 'FontWeight', 'bold', 'Interpreter', 'latex');
    axis square;                             % square plot

    hold off;
end  

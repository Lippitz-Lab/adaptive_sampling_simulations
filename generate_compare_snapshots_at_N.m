% Generate comparable residual snapshots for square, Daemon, and adaptive sampling.
% Match target MSE values using reproducible noise realizations.

function generate_compare_snapshots_at_N(idmol)

    % Preserve the caller's random-number-generator state.
    previousRng = rng;
    cleanupRng = onCleanup(@() rng(previousRng));  

    % Load the molecule data and previously calculated scan statistics.
    load("molecules.mat", "trace_r", "gamma0", "t1", "t3");
    load("algo_chi.mat", "allchis");
    load("daemon_and_square.mat", "noiselevels");

    % Use a reproducible random sequence for this calculation.
    rng(1);



    % Maximum number of measurements in the scan simulations.
    Nmax = 1e3;

    % Target normalized mean-square errors for the three methods.
    N_show = [1e3, 1e4, 1e5];
    nTargets = numel(N_show);

    % Use the precomputed curves to estimate suitable noise levels.
    % Rows correspond to the daemon and adaptive algorithms.
    allchis_selected = [
        allchis.daemon(:).';
        allchis.algo(:).'
    ];


    % Extract the selected molecule and its total signal power.
    signal_truth = squeeze(trace_r(idmol, :, :));
    signalNorm = sum(abs(signal_truth).^2, "all");

    % Allocate one panel for each method and target MSE.
    panels = repmat(struct( ...
        "coord", [], ...
        "window", [], ...
        "Nmeas", [], ...
        "N", [], ...
        "chitarget", []), 1, 3 * nTargets);

    panelId = 0;

    for id = 1:nTargets

        % -------------------------------------------------------------
        % Square scan
        % -------------------------------------------------------------

        % The square method measures a fixed 7-by-7 region.
        edgelength = 7;
        window = zeros(size(signal_truth));
        window(1:edgelength, 1:edgelength) = 1;

        Nmeas = N_show(id) ./ sum(window(:)) .* window;

        panelId = panelId + 1;
        panels(panelId).coord = [1, id];
        panels(panelId).window = window;
        panels(panelId).Nmeas = Nmeas;
        panels(panelId).N = sum(Nmeas, "all");
        panels(panelId).chitarget = NaN;

        % -------------------------------------------------------------
        % Daemon scan
        % -------------------------------------------------------------


        Nmax = 1000;
        noiselevel = 1 ./ sqrt(N_show(id) / 1000);

        [~, ~, snapshot] = daemon_scan( ...
                signal_truth, noiselevel, Nmax, 0);

        panelId = panelId + 1;
        panels(panelId).coord = [2, id];
        panels(panelId).window = snapshot.window;
        panels(panelId).Nmeas = snapshot.Nmeas ./ noiselevel.^2;
        panels(panelId).N = sum(panels(panelId).Nmeas, "all");
        panels(panelId).chitarget = NaN;

        % -------------------------------------------------------------
        % Adaptive algorithm
        % -------------------------------------------------------------


        Nmax = 1000;
        noiselevel = 1 ./ sqrt(N_show(id) / 1000);

        [~, ~, ~, snapshot] = adaptive_algo_scan( ...
                signal_truth, noiselevel, Nmax, gamma0, ...
                t1 + t3, 0);


        panelId = panelId + 1;
        panels(panelId).coord = [3, id];
        panels(panelId).window = snapshot.window;
        panels(panelId).Nmeas = snapshot.Nmeas ./ noiselevel.^2;
        panels(panelId).N = sum(panels(panelId).Nmeas, "all");
        panels(panelId).chitarget = NaN;
    end


    % create noise realizations for display
    rng(1);
    noise_r0 = randn(size(signal_truth)) + ...
        1i .* randn(size(signal_truth))./ sqrt(2);

    % Generate an independent noise realization for non-rephasing data.
    noise_nr0 = randn(size(signal_truth)) + ...
        1i .* randn(size(signal_truth))./ sqrt(2);

    % Save the sampling panels and the selected noise realizations.
    save("compare_snapshots_at_N.mat", ...
        "panels", "noise_r0", "noise_nr0", "N_show");

end
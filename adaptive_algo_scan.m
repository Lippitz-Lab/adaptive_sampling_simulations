% simulates a scan with the adaptive algorithm
% We use random numbers to simulate shot noise

function [chilist, mlist, param_trace, snapshot] = adaptive_algo_scan( ...
        signal_truth, noiselevel, Nmax, gamma0, times, chi_snapshot)

    % Keep the full complex value
    signal_truth = squeeze(signal_truth);

    % Initialize
    signalNorm = sum(abs(signal_truth).^2, "all");
    noiseVariance = noiselevel^2;

    snapshot = [];
    dataSum = zeros(size(signal_truth), "like", signal_truth);

    param_trace = NaN(Nmax, 4);
    chilist = NaN(Nmax, 1);
    mlist = NaN(Nmax, 1);

    % Three starting points
    Nmeas = zeros(size(signal_truth));
    Nmeas(1:2, 1) = 1;
    Nmeas(1, 2) = 1;

    measuredPoints = Nmeas > 0;

    % First three measurements by hand
    initialData = signal_truth + ...
        (randn(size(signal_truth)) + 1i .* randn(size(signal_truth))) ...
        .* noiselevel / sqrt(2);

    dataSum(measuredPoints) = initialData(measuredPoints);
    noisestore = NaN(Nmax, 1);

    for run = 4:Nmax

        % Current mean estimate from measured data
        data = NaN(size(signal_truth), "like", signal_truth);
        data(measuredPoints) = ...
            dataSum(measuredPoints) ./ Nmeas(measuredPoints);

        % Fit an exponential decay to measured |c'_n|^2
        [cn2, p] = fit_exp_decay(abs(data).^2, times, Nmeas, gamma0, noiselevel);

        param_trace(run, :) = p;

        % Calculate window function based on surrogate model
        window = cn2 .* Nmeas ./ ...
            (cn2 .* Nmeas + noiseVariance);

        % Determine chi-square from expected residuals
        residuals = abs(signal_truth).^2;
        residuals(measuredPoints) = ...
            window(measuredPoints) .* noiseVariance ./ ...
            Nmeas(measuredPoints);

        chisq = sum(residuals, "all") ./ signalNorm;

        chilist(run) = chisq;
        mlist(run) = sum(measuredPoints, "all");

        % Store snapshot once, as soon as chi_snapshot is passed
        if isempty(snapshot) && chisq < chi_snapshot
            snapshot.Nmeas = Nmeas;
            snapshot.window = window;
            snapshot.residuals_sq = residuals ./ signalNorm;
            snapshot.chisq = chisq;
            snapshot.N = sum(Nmeas, "all");
        end

        % Determine next point to measure
        delta = cn2 ./ (cn2 .* Nmeas + noiseVariance) - ...
            cn2 ./ (cn2 .* (Nmeas + 1) + noiseVariance);

        [~, idnext] = max(delta, [], "all");

        % Measure the selected point
        randomness =  (randn + 1i .* randn) .* noiselevel ./ sqrt(2);
        measurement = signal_truth(idnext) + randomness;
        noisestore(run) = abs(randomness).^2;


        Nmeas(idnext) = Nmeas(idnext) + 1;
        dataSum(idnext) = dataSum(idnext) + measurement;
        measuredPoints(idnext) = true;
    end

    % store snapshot at the end if chi_snapshot == 0
    if  chi_snapshot == 0
        snapshot.Nmeas = Nmeas;
        snapshot.window = window;
        snapshot.residuals_sq = residuals ./ signalNorm;
        snapshot.chisq = chisq;
        snapshot.N = sum(Nmeas, "all");
    end


end
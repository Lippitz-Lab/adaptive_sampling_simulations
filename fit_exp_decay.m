% Fit an exponentially decaying function to 2D data.
% The decay rate and amplitude are constrained to a physically meaningful range.

function [model, param] = fit_exp_decay(data0, time0, Nmeas, gamma0, noiselevel)


    % we roughly know the amplitude of the first coefficient
    % so we shift the time axis so that t=0 at this coeff
    time0 = time0 - min(time0(:));

    % Fit only measured points
    mask = Nmeas > 0;
    data = data0(mask);
    time = time0(mask);
    Nmeas = Nmeas(mask);


    noiseVariance =  noiselevel^2 ./ Nmeas;

    [~, id] = min(time(:));
    amp00 = max( [data(id) - noiseVariance(id), 20 .* realmin]);


    % A rough fit is sufficient here. The amplitude grid is centered on the
    % expected value of one and includes the physically allowed lower bound.
   % bound_factor = 2;
    bound_interval = [0.5, 2 ];
    gammaBounds = log10(bound_interval .* gamma0);
    gammaGrid = logspace(gammaBounds(1), gammaBounds(2), 11);

    ampBound = log10(bound_interval .* amp00);
    amplitudeGrid = logspace(ampBound(1), ampBound(2), 11);

    % Arrange dimensions as amplitude x gamma x measurement. I
    time = reshape(time, 1, 1, []);
    noiseVariance = reshape(noiseVariance, 1, 1, []);
    scaledData = reshape(data, 1, 1, []) ./ noiseVariance;
    basis = exp(-2 .* time .* reshape(gammaGrid, 1, [], 1));

    % The noncentrality parameter is the expected squared signal divided by
    % the noise variance.
    noncentrality = reshape(amplitudeGrid, [], 1, 1) .* basis ...
                  ./ noiseVariance;

    % This matches the manuscript equation for the 2D noncentral chi-square
    % likelihood:
    %
    %   L(A,gamma) = sum_n [ s_n^2 / (sigma^2 / N_n)
    %                      - log(I0( 2 * sqrt(s_n^2 |X_n|^2)
    %                                 / (sigma^2 / N_n) )) ]
    %
    argument = 2 .* sqrt(noncentrality .* scaledData);

    % The scaled Bessel function besseli(..., 1) avoids overflow; adding the
    % argument recovers log(I0(z)) for z >= 0.
    logBessel = log(max(besseli(0, argument, 1), realmin)) + argument;

    logLikelihood = sum(logBessel - noncentrality, 3);

    % We  maximize the log-likelihood. 
    [~, bestIndex] = max(logLikelihood(:));
    [bestAmplitudeIndex, bestGammaIndex] = ind2sub( ...
        size(logLikelihood), bestIndex);

    amp =  amplitudeGrid(bestAmplitudeIndex);
    gamma =  gammaGrid(bestGammaIndex);

    % Calculate model over the complete time range
    model =   amp .* exp(-time0 .* 2.* gamma);

    [~, id] = min(time0(:));

    param = [amp ./ amp00, gamma, model(id) + noiseVariance(id), data0(id)];

end

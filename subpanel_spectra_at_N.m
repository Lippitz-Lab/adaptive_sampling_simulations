% Compare measured and reference 2D spectra and their diagonal cuts.
% Display the spectra, calculate normalized chi-square errors, and plot residual differences.

function [diagcut, diagChisq] = subpanel_spectra_at_N( ...
        t_r0, t_nr0, Nmeas, window, noise_r, noise_nr, ...
        t_list, col, row, layout)

    % Plot and analysis settings.
    plotCfg.fftPadding = 10;
    plotCfg.contourLevels = -10:10;
    plotCfg.spectrumClim = [-9, 0.4];
  

    cmap =brewermapNL(500, '-RdBu',1.2);
    mygray = [1,1,1] .* 0.8;
    myblue = cmap(50,:);
    myred = cmap(430,:);

    % Simulate the measured rephasing trace only at measured points.
    measuredPoints = Nmeas > 0;
    t_r = zeros(size(t_r0), "like", t_r0);
    t_r(measuredPoints) = window(measuredPoints) .* ...
        (t_r0(measuredPoints) + ...
        noise_r(measuredPoints) ./ sqrt(Nmeas(measuredPoints)));

    % Simulate the measured non-rephasing trace only at measured points.
    t_nr = zeros(size(t_nr0), "like", t_nr0);
    t_nr(measuredPoints) = window(measuredPoints) .* ...
        (t_nr0(measuredPoints) + ...
        noise_nr(measuredPoints) ./ sqrt(Nmeas(measuredPoints)));

    % Fourier-transform measured and reference traces.
    [spec_r, ~] = my2Dfft_2axes_Trapezoidal( ...
        t_r, t_list, t_list, plotCfg.fftPadding);
    [spec_r0, ~] = my2Dfft_2axes_Trapezoidal( ...
        t_r0, t_list, t_list, plotCfg.fftPadding);

    [spec_nr, ~] = my2Dfft_2axes_Trapezoidal( ...
        t_nr, t_list, t_list, plotCfg.fftPadding);
    [spec_nr0, omega_fft] = my2Dfft_2axes_Trapezoidal( ...
        t_nr0, t_list, t_list, plotCfg.fftPadding);

    % Apply rotating-frame frequency offsets.
    % The current zero offsets retain the FFT frequency coordinates.
    w1t = 0 + omega_fft;
    w3t = 0 - omega_fft;

    % Construct purely absorptive measured and reference spectra.
    pure_abs = real(spec_r + flipud(spec_nr));
    pure_abs0 = real(spec_r0 + flipud(spec_nr0));

    % Plot the 2D spectrum.
    axSpectrum = nexttile(tilenum(layout, row, col));

  %  imagesc(axSpectrum, w3t, w1t, pure_abs);
  %  hold(axSpectrum, "on");

    contourf(axSpectrum, w3t, w1t, pure_abs, ...
        plotCfg.contourLevels, ...
            "LineColor", 'none');

    
    %  "LineColor", mygray,'LineWidth',0);

  %  hold(axSpectrum, "off");
    set(axSpectrum, "YDir", "normal");
    axis(axSpectrum, "square");
    set(axSpectrum, "XTickLabel", [], "YTickLabel", []);

    colormap(axSpectrum, cmap);
    clim(axSpectrum, [plotCfg.spectrumClim(1), -plotCfg.spectrumClim(1)] )

    spectrumChisq = sum(abs(pure_abs - pure_abs0).^2, "all") ./ ...
            sum(abs(pure_abs0).^2, "all");
  
    text(axSpectrum, 1, 0.9, ...
        sprintf("\\chi^2=%.1e ", spectrumChisq), ...
        "Color", "k", ...
        "HorizontalAlignment", "right", ...
        "Units", "normalized");

    % Extract diagonal cuts from measured and reference spectra.
    diagcut = diag(fliplr(pure_abs));
    diagcut0 = diag(fliplr(pure_abs0));


   diagChisq = sum(abs(diagcut - diagcut0).^2, "all") ./ ...
            sum(abs(diagcut0).^2, "all");
  
    % Plot the region between the measured and reference diagonal cuts.
    axDiagonal = nexttile(tilenum(layout, row, col + 2));

    x = 1:numel(diagcut);
    plot(x,x .* 0,'k:')

    hold on
    fill(axDiagonal, ...
        [x, fliplr(x)], ...
        [diagcut(:).', fliplr(diagcut0(:).')], ...
        myred, ...
        "LineStyle", "none");
 

    xlim(axDiagonal, [1, numel(diagcut)]);
    ylim(axDiagonal, plotCfg.spectrumClim);
    set(axDiagonal, "XTickLabel", [], "YTickLabel", []);
    axis(axDiagonal, "square");

    text(axDiagonal, 1, 0.1, ...
        sprintf("\\chi^2=%.1e ", diagChisq), ...
        "Color", "k", ...
        "HorizontalAlignment", "right", ...
        "Units", "normalized");
 
end
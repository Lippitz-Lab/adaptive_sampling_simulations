% Fig 2 in the paper
% plot linear spectrum, time domain data and cuts, and 2d spectrum

function fig_mol_spectra(idmol, use_export_fig)
    

    % Load simulated data
    load("molecules.mat", ...
         "Nmol", ...
        "spec_lin", "omega_lin", ...
        "trace_r", "trace_nr", ...
        "t1", "t3", "t_list", ...
        "gamma0", "omega21", "omega43");

    % Plot configuration
    plotCfg.pmax = 6;
    plotCfg.linearXLim = [2.12, 2.48];
    plotCfg.timeXLim = [0, 120];
    plotCfg.timeYLim = [0, 120];
    plotCfg.lineCutXLim = [0, 200];
    plotCfg.spectrumClim = [-9, 0.4];
    plotCfg.contourLevels = -10:0;
    plotCfg.fftpadding = 10;

    % demo molecule
    demo = squeeze(trace_r(idmol, :, :));
    demo_nr = squeeze(trace_nr(idmol, :, :));


    cmap =brewermapNL(500, '-RdBu',1.2);
    mygray = [1,1,1] .* 0.8;
    myblue = cmap(50,:);
    myred = cmap(430,:);

    fig = figure;
    tiledlayout(2, 2,'Padding','tight')

    % Linear spectra
    ax1 = nexttile;

    ensemble = squeeze(mean(spec_lin, 1));

    area(omega_lin, ensemble,'FaceColor',mygray,'EdgeColor','none')
    hold on
    plot(omega_lin, spec_lin(idmol, :),'LineWidth',1.5,'Color',myred)
    plot(omega_lin, spec_lin(idmol+1, :),'Color',myblue)
    plot(omega_lin, spec_lin(idmol+2, :),'Color',myblue)
    plot(omega_lin, spec_lin(idmol+3, :),'Color',myblue)

    xticks([2.2, 2.3, 2.4])
   

    xlabel("frequency (rad/fs)")
    ylabel("linear absorption")
    axis square
    xlim(plotCfg.linearXLim)
    yticks(0)

    % Time domain
    ax2 = nexttile;

    pmax = plotCfg.pmax;
    traceImage = abs(demo(1:pmax, 1:pmax));

    scaleddata = log10(max(traceImage, eps));
    scaleddata = scaleddata - min(scaleddata(:));

    imagesc( ...
        t_list(1:pmax), ...
        t_list(1:pmax), ...
        scaleddata)

    peakval = max(scaleddata(:));

    set(gca, "YDir", "normal")
    axis square
    
    colormap(ax2, cmap);
    %colormap summer
    clim([-peakval, peakval] )

    cbh = colorbar;
    cbh.Ticks = peakval -[2,1,0];
    cbh.Limits = [0, peakval + 0.1];
    cbh.TickLabels = ["10^{-2}", "10^{-1}", "1"];
    cbh.Label.String = "rephasing |c_n|";

    xticks([0, 50, 100]);
    yticks([0, 50, 100]);
    xlim(plotCfg.timeXLim)
    ylim(plotCfg.timeYLim)
    xlabel("delay t_{3} (fs)")
    ylabel("delay t_{1} (fs)")

    % Single rephasing line cuts
    ax3 = nexttile;

    time = t1 + t3;

    diagtime = diag(time);

    stairs(time(1,:) - time(1,1)/2, abs(demo( 1, :)), ...
        "LineWidth", 1, "Color", myblue)
    hold on
  %  plot(diag(time), abs(diag(demo)),'o')
    stairs(diagtime - diagtime(1), abs(diag(demo)), ...
        "LineWidth", 1, "Color", myred)

    plot(t_list, 2 .* exp(-t_list .* gamma0), "--", ...
        "Color", "k")

    xlabel("delay t_1 + t_3 (fs)")
    ylabel("rephasing |c_n|")
    xlim(plotCfg.lineCutXLim)
    axis square
    set(gca,'YScale','log');

    % Single PA spectrum
   ax4 =  nexttile;

    [spec_r, ~] = my2Dfft_2axes_Trapezoidal( ...
        demo, t_list, t_list, plotCfg.fftpadding);

    [spec_nr, omega_fft] = my2Dfft_2axes_Trapezoidal( ...
        demo_nr, t_list, t_list, plotCfg.fftpadding);

    w1t = omega21 + omega_fft;
    w3t = omega43 - omega_fft;

    pure_abs = real(spec_r + flipud(spec_nr));

    % imagesc(w3t, w1t, pure_abs)
    % hold on
    % contour(w3t, w1t, pure_abs, plotCfg.contourLevels, ...
    %     "LineColor", mygray)
    % hold off
    contourf(ax4, w3t, w1t, pure_abs, plotCfg.contourLevels, ...
        "LineColor", "none")

    colormap(ax4, cmap);
    clim([plotCfg.spectrumClim(1), -plotCfg.spectrumClim(1)] )


    set(gca, "YDir", "normal")
    axis square
    cbh = colorbar;
    cbh.Label.String = "real purely absorptive";
    cbh.Limits = plotCfg.spectrumClim;

    xticks([2.2, 2.3, 2.4])
    yticks([2.2, 2.3, 2.4])
    xlabel("frequency \omega_{3} (rad/fs)")
    ylabel("frequency \omega_{1} (rad/fs)")

    scifig( ...
        "width", 8.2, ...
        "height", 6, ...
        "font", "Arial", ...
        "fontsize", 7, ...
        "tickdir", "in", ...
        "linewidth", 0.5);
    % Set the target dimensions on the figure object
   % fig.Units = "centimeters";
   % fig.Position(3:4) = [8.2, 5]; % [width, height]
 

    if use_export_fig
        export_fig("fig_mol_spectra.pdf",'-pdf','-transparent','-r300','-nocrop','-silent');
    else
        exportgraphics(fig, "fig_mol_spectra.pdf","ContentType", "vector"); 
    end


end
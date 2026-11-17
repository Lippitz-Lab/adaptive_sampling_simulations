% Figure chi over N, Fig. 3 in the paper
%
function fig_chi_over_N

    cmap =brewermapNL(500, '-RdBu',1.2);
    mygray = [1,1,1] .* 0.8;
    myblue = cmap(50,:);
    myred = cmap(430,:);
    mygreen = [0, 160, 62] ./ 255 .* 1.5;




    load("algo_chi.mat", "allchis");

    % Plot configuration
    plotCfg.gray = 0.5 .* [1, 1, 1];
    plotCfg.lineWidth = 0.5;
    plotCfg.mainLineWidth = 1;
    plotCfg.xLim = [1e1, 1e8];
    plotCfg.yLim = [1e-4, 7];
    plotCfg.xTicks = 10.^(1:8);
    plotCfg.figureWidth = 8.2;
    plotCfg.figureHeight = 5;

    
    x = allchis.x(:);
    algoMean = allchis.algo(:);
    daemonMean = allchis.daemon(:);

    shrinkfactor = 0.25;
    algoUpper = algoMean + shrinkfactor * allchis.algo_std(:);
    algoLower = max(algoMean - shrinkfactor* allchis.algo_std(:), eps);

    daemonUpper = daemonMean + shrinkfactor* allchis.daemon_std(:);
    daemonLower = max(daemonMean - shrinkfactor* allchis.daemon_std(:), eps);

   load("daemon_and_square.mat",  "allchis");

    sq3Upper = allchis.square3(:) + shrinkfactor* allchis.square3_std(:);
    sq3Lower = allchis.square3(:) - shrinkfactor* allchis.square3_std(:);

   sq5Upper = allchis.square5(:) + shrinkfactor* allchis.square5_std(:);
    sq5Lower = allchis.square5(:) - shrinkfactor* allchis.square5_std(:);

     sq7Upper = allchis.square7(:) + shrinkfactor* allchis.square7_std(:);
    sq7Lower = allchis.square7(:) - shrinkfactor* allchis.square7_std(:);

    fig = figure;

    loglog(x, algoMean, ...
        "LineWidth", plotCfg.mainLineWidth, "Color", myred);
    
    hold on

    loglog(x, allchis.square5(:), ...
        "LineWidth", plotCfg.lineWidth, "Color", "k");
    loglog(x, allchis.square3(:), ...
        "LineWidth", plotCfg.lineWidth, "Color", "k");
    loglog(x, allchis.square7(:), ...
        "LineWidth", plotCfg.lineWidth, "Color", "k");

    loglog(x, daemonMean, ...
        "LineWidth", plotCfg.mainLineWidth, "Color", myblue);

    %      Transparent uncertainty bands
    fill([x; flipud(x)], [algoUpper; flipud(algoLower)], ...
        'r', "FaceAlpha", 0.15, "EdgeColor", "none");

    fill([x; flipud(x)], [daemonUpper; flipud(daemonLower)], ...
        'b', "FaceAlpha", 0.15, "EdgeColor", "none");

    fill([x; flipud(x)], [sq3Upper; flipud(sq3Lower)], ...
        'k', "FaceAlpha", 0.15, "EdgeColor", "none");

    fill([x; flipud(x)], [sq5Upper; flipud(sq5Lower)], ...
        'k', "FaceAlpha", 0.15, "EdgeColor", "none");

   fill([x; flipud(x)], [sq7Upper; flipud(sq7Lower)], ...
        'k', "FaceAlpha", 0.15, "EdgeColor", "none");

        xlim(plotCfg.xLim);
    ylim(plotCfg.yLim);
    xticks(plotCfg.xTicks);

    text(2e5, 5e-4, "Daemon", ...
        "FontSize", 5, "FontName", "Arial");
    text(2e6, 1.5e-1, "fixed square", ...
        "FontSize", 5, "FontName", "Arial");

    xlabel("total integration time T / T_0");
    ylabel("normalized MSE");

    hold off;

    scifig( ...
        "width", plotCfg.figureWidth, ...
        "height", plotCfg.figureHeight, ...
        "font", "Arial", ...
        "fontsize", 7, ...
        "tickdir", "in", ...
        "linewidth", plotCfg.lineWidth);
    
    % Set dimensions directly on the figure object
    fig.Units = "centimeters";
    fig.Position(3:4) = [plotCfg.figureWidth, plotCfg.figureHeight];
    
    % Export without size options
    exportgraphics(fig, "fig_chi_over_N.pdf", ...
        "ContentType", "vector");

    exportgraphics(fig, "fig_chi_over_N.png")

end
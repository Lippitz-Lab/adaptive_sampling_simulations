% residuals subpanel in Fig .4

function maxres = subpanel_residuals_at_N(panel, layout)

    plotSize = 10;

    ax = nexttile(tilenum(layout, panel.coord(2), panel.coord(1)));

    % Identify measured pixels and calculate measurement counts.
    touched = panel.Nmeas > 0;
    m = sum(touched, "all");
    N = sum(panel.Nmeas, "all");

    % Calculate the outline of the measured region.
    myborder = bwperimtrace(touched);

    residualImage = abs(panel.res(1:plotSize, 1:plotSize));
    imagesc(residualImage);
    hold on

    % Overlay the measured-region boundary when available.
    if ~isempty(myborder) 
        for id =1:numel(myborder)
           plot(myborder{id}(:, 1), myborder{id}(:, 2), ...
                "Color", 'k', "LineWidth", 1,'LineStyle',':');
        end
    end

    set(gca, "YDir", "normal");
    axis square
    set(gca, "XTickLabel", []);
    set(gca, "YTickLabel", []);

    % Add summary information to the upper-right corner.
   text(plotSize, plotSize, ...
        sprintf("\\chi^2 = %.1e", panel.realchisq), ...
        "Color", "k", "HorizontalAlignment", "right");

    % text(plotSize, plotSize -1, ...
    %     "T / T_0 = " + commaformat(round(N)), ...
    %     "Color", "w", "HorizontalAlignment", "right");

    % text(plotSize, plotSize - 2, ...
    %     "m = " + commaformat(m), ...
    %     "Color", "w", "HorizontalAlignment", "right");


    hold off

    maxres = max(residualImage, [], "all");

end


% Format a nonnegative integer with thousands separators.
function text = commaformat(value)
   text = string(sprintf("%.0f", value));
   text = regexprep(text, "(\d)(?=(\d{3})+$)", "$1 ");
end
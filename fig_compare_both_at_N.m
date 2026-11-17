% Figure 4 in the paper
% Compare residuals from square, daemon, and adaptive sampling strategies.

function fig_compare_both_at_N(idmol, use_export_fig)

    cmap =brewermapNL(500, '-RdBu',1.2);
    mygray = [1,1,1] .* 0.8;
    myblue = cmap(50,:);
    myred = cmap(430,:);



    % Load the simulated signal and the selected molecule index.
    load("molecules.mat", "trace_r", "trace_nr", "t_list");

    % Load sampling patterns and the noise realization used for comparison.
    load("compare_snapshots_at_N.mat", ...
        "panels", "noise_r0", "noise_nr0", "N_show");

    % Extract the selected molecule's rephasing signal.
    signal_truth = squeeze(trace_r(idmol, :, :));
    signalNorm = sum(abs(signal_truth).^2, "all");

    % Create a 3-by-7 layout:
    % rows correspond to target chi-square values;
    % columns correspond to the three sampling methods + spectra
    fig = figure;
    layout = tiledlayout(fig, 3, 7,'TileSpacing','compact','Padding','tight');

    % Store the maximum residual value for each panel.
    nPanels = numel(panels);
    maxres = NaN(1, nPanels);

    for id = 1:nPanels

        % Add noise only at positions that were measured.
        measured = zeros(size(signal_truth), "like", signal_truth);
        measuredPoints = panels(id).Nmeas > 0;

        measured(measuredPoints) = ...
            panels(id).window(measuredPoints) .* ...
            noise_r0(measuredPoints) ./ ...
            sqrt(panels(id).Nmeas(measuredPoints));

        % Combine the unmeasured signal and measured noisy data.
        panels(id).res = ...
            (1 - panels(id).window) .* signal_truth + measured;

        % Calculate the normalized residual intensity.
        panels(id).realchisq = ...
            sum(abs(panels(id).res).^2, "all") ./ signalNorm;

        % Draw the residual map and record its maximum value.
        maxres(id) = subpanel_residuals_at_N(panels(id), layout);
    end

    % Sampling methods shown in the columns.
    mytitles = ["square7", "daemon", "algo"];

    for col = 1:3
        nexttile(tilenum(layout, 1, col));
        title(mytitles(col));
    end

    % Group the panels by chitarget and find one color limit per chitarget.
    maxres = max(reshape(maxres, 3, 3), [], 1);

    for id = 1:numel(N_show)

        % Label each row with its target chi-square value.
        nexttile(tilenum(layout, id, 1));
        ylabel(sprintf("T = %s T_0", commaformat(N_show(id))));

        % Use consistent color limits across each row.
        for col = 1:3
            ax = nexttile(tilenum(layout, id, col));
            colormap(ax,cmap)
            clim([-maxres(id), maxres(id)]);
        end
    end

%-----------
  % Extract the selected molecule without overwriting the loaded arrays.
    traceR = squeeze(trace_r(idmol, :, :));
    traceNR = squeeze(trace_nr(idmol, :, :));

    % Select the square7 and adaptive panels for each target chi-square.
    idlist = [(1:3:7), (3:3:9)];

  

    for idt = 1:numel(idlist)

        id = idlist(idt);

        % Map adaptive panels from the third method column to the
        % second displayed column. The panel ordering remains unchanged.
        x = panels(id).coord(1);
        if x == 3
            x = 2;
        end

        subpanel_spectra_at_N( ...
            traceR, traceNR, ...
            panels(id).Nmeas, ...
            panels(id).window, ...
            noise_r0, noise_nr0, ...
            t_list, ...
            x+3, ...
            panels(id).coord(2), ...
            layout);
    end


       % Label the displayed method columns.
    mytitles = ["square7", "algo", "square7", "algo"];

    for col = 1:4
        nexttile(tilenum(layout, 1, col+3));
        title(mytitles(col));
    end

%------------

  scifig( ...
        "width", 8.2 * 2, ...
        "height", 7, ...
        "font", "Arial", ...
        "fontsize", 5, ...
        "tickdir", "in", ...
        "linewidth", 0.5);

    if use_export_fig
       export_fig(fig, "fig_compare_both_at_N",'-pdf','-nocrop','-silent');
       export_fig(fig, "fig_compare_both_at_N",'-svg','-painters','-nocrop','-silent');
    else
         exportgraphics(fig, "fig_compare_both_at_N.pdf", ...
            "ContentType", "vector"); 
    end


end

% Format a nonnegative integer with thousands separators.
function text = commaformat(value)
   text = string(sprintf("%.0f", value));
   text = regexprep(text, "(\d)(?=(\d{3})+$)", "$1 ");
end
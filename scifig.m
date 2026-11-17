function scifig(options)
% Format figures for scientific publications.
% 
%
% based on code by
% Urban J. Wuensch - wuensch@chalmers.se
% Chalmers University of Technology
% Department of Architecture and Civil Engineering
% Sven Hultins Gata 6
% 41296 Gothenburg (Sweden)


arguments
    options.width (1,1) {mustBeNumericOrLogical, mustBePositive} = 8
    options.height (1,1) {mustBeNumericOrLogical, mustBePositive} = 7
    options.font {mustBeText} = "Arial"
    options.fontsize (1,1) {mustBeNumericOrLogical, mustBePositive} = 9
    options.axes {mustBeText, ...
        mustBeMember(options.axes, ["create", "nocreate"])} = "nocreate"
    options.tickdir {mustBeText, ...
        mustBeMember(options.tickdir, ["both", "in", "out", "none"])} = "both"
    options.linewidth (1,1) {mustBeNonnegative} = 0.5
end

% Use the current figure, preserving its existing position.
fig = gcf;

% Configure figure appearance and size in centimeters.
fig.Units = "centimeters";
fig.InvertHardcopy = "off";
fig.Color = [1, 1, 1];

position = fig.Position;
fig.Position = [position(1:2), options.width, options.height];

movegui(fig, "center");

% Find existing axes without creating new ones.
axesHandles = findobj(fig, "Type", "axes");

if isempty(axesHandles)
    if options.axes == "create"
        axesHandles = axes("Parent", fig);
    else
        % Apply figure-wide formatting even if no axes exist.
        fontsize(fig, options.fontsize, "points");
        fontname(fig, options.font);
        return
    end
end

% Apply figure-wide font settings in one call each.
fontsize(fig, options.fontsize, "points");
fontname(fig, options.font);

% Apply common properties to all axes at once.
set(axesHandles, ...
    "TickDir", options.tickdir, ...
    "LineWidth", options.linewidth);

% Tiled layouts manage their own outer positions.
for ax = axesHandles(:).'
    if ~strcmp(ax.Parent.Type, "tiledlayout")
        outerPosition = ax.OuterPosition;
        outerPosition(outerPosition < 0) = 0;
        ax.OuterPosition = outerPosition;
    end
end

drawnow

end
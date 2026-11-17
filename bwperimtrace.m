function cellout = bwperimtrace(in,xlims,ylims)
%
% based on code by
% Daniel Warren 10/11/2014, Department of Oncology, University of Oxford

% Trace the outer boundaries of nonzero pixels using 4-connectivity.
%
% Coordinates describe pixel edges. The optional xlims and ylims specify
% the physical locations of the centers of the first and last pixels.
%
% Note: this function transposes the input internally. Call it with the
% original image orientation; do not transpose the input at the call site.

    arguments
        in {mustBeNumericOrLogical, mustBeNonempty, mustBeFinite}
        xlims (1,2) double {mustBeFinite} = [1, size(in, 1)]
        ylims (1,2) double {mustBeFinite} = [1, size(in, 2)]
    end

    if ~ismatrix(in)
        error("in must be a nonempty two-dimensional array.");
    end

    if any(xlims(1) == xlims(2)) && size(in, 1) > 1
        error("xlims must contain two distinct values.");
    end

    if any(ylims(1) == ylims(2)) && size(in, 2) > 1
        error("ylims must contain two distinct values.");
    end

    % Convert numeric input to a logical image.
    in = in ~= 0;

    % Return an empty cell array when no pixels are active.
    if ~any(in(:))
        cellout = cell(0, 1);
        return
    end

    % The implementation uses x as the first matrix dimension.
    in = in.';

    % Use linspace so the coordinate vectors have exactly the required
    % number of entries, including for one-row or one-column images.
    if size(in, 1) == 1
        xpix = xlims(1);
        dx = 1;
    else
        xpix = linspace(xlims(1), xlims(2), size(in, 1));
        dx = (xlims(2) - xlims(1)) / (size(in, 1) - 1);
    end

    if size(in, 2) == 1
        ypix = ylims(1);
        dy = 1;
    else
        ypix = linspace(ylims(1), ylims(2), size(in, 2));
        dy = (ylims(2) - ylims(1)) / (size(in, 2) - 1);
    end

    % Pixel boundaries are half a pixel away from pixel centers.
    xbnd = [xpix - 0.5 * dx, xpix(end) + 0.5 * dx].';
    ybnd = [ypix - 0.5 * dy, ypix(end) + 0.5 * dy].';

    [indx, indy] = find(in);
    nPixels = numel(indx);

    % Store four boundary segments per active pixel.
    points = zeros(8 * nPixels, 2);

    x1 = xbnd(indx) ;
    x2 = xbnd(indx + 1);
    y1 = ybnd(indy);
    y2 = ybnd(indy + 1);

    points(1:8:end, :) = [x1, y1];
    points(2:8:end, :) = [x2, y1];
    points(3:8:end, :) = points(2:8:end, :);
    points(4:8:end, :) = [x2, y2];
    points(5:8:end, :) = points(4:8:end, :);
    points(6:8:end, :) = [x1, y2];
    points(7:8:end, :) = points(6:8:end, :);
    points(8:8:end, :) = points(1:8:end, :);

    % Remove pairs of opposing segments shared by neighboring pixels.
    midpts = 0.5 .* (points(1:2:end, :) + points(2:2:end, :));
    dirs = points(1:2:end, :) - points(2:2:end, :);

    removeSegment = false(size(midpts, 1), 1);

    upinds = find(dirs(:, 2) > 0);
    downinds = find(dirs(:, 2) < 0);
    [~, upRemove, downRemove] = intersect( ...
        midpts(upinds, :), midpts(downinds, :), "rows");

    removeSegment(upinds(upRemove)) = true;
    removeSegment(downinds(downRemove)) = true;

    rightinds = find(dirs(:, 1) > 0);
    leftinds = find(dirs(:, 1) < 0);
    [~, leftRemove, rightRemove] = intersect( ...
        midpts(leftinds, :), midpts(rightinds, :), "rows");

    removeSegment(leftinds(leftRemove)) = true;
    removeSegment(rightinds(rightRemove)) = true;

    % Each segment consists of two consecutive rows.
    keepPoint = repelem(~removeSegment, 2);
    points = points(keepPoint, :);
    
    % Reorder connected segments into closed boundary traces.
    
    nSegments = size(points, 1) / 2;
    segmentStart = points(1:2:end, :);
    segmentEnd = points(2:2:end, :);
    unused = true(nSegments, 1);
    
    cellout = cell(nSegments, 1);
    nObjects = 0;
    
    while any(unused)
    
        % Start a new boundary with the first unused segment.
        segmentId = find(unused, 1, "first");
        startPoint = segmentStart(segmentId, :);
        currentPoint = segmentEnd(segmentId, :);
    
        boundary = [
            segmentStart(segmentId, :)
            segmentEnd(segmentId, :)
        ];
    
        unused(segmentId) = false;
        nUsed = 1;
    
        while ~isequal(currentPoint, startPoint)
    
            % Prevent an endless loop if the boundary cannot be closed.
            if nUsed >= nSegments
                error("Boundary tracing failed: boundary is not closed.");
            end
    
            % Prefer a segment whose start matches the current endpoint.
            candidate = find( ...
                unused & all(segmentStart == currentPoint, 2), ...
                1, "first");
    
            if ~isempty(candidate)
                nextPoint = segmentEnd(candidate, :);
                boundary(end + 1, :) = nextPoint;
            else
                % Otherwise, reverse a segment whose end matches.
                candidate = find( ...
                    unused & all(segmentEnd == currentPoint, 2), ...
                    1, "first");
    
                if isempty(candidate)
                    error("Boundary tracing failed: disconnected segment.");
                end
    
                nextPoint = segmentStart(candidate, :);
                boundary(end + 1, :) = nextPoint;
            end
    
            unused(candidate) = false;
            currentPoint = nextPoint;
            nUsed = nUsed + 1;
        end
    
        nObjects = nObjects + 1;
        cellout{nObjects} = boundary;
    end
    
    cellout = cellout(1:nObjects);

end
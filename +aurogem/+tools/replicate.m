% Description:
%   aaa
%
% Example usage:
%   [bc, resnorm, E2_bg, E3_bg, v2_int, v3_int, weight0, bound] = replicate(track, image, xg, opts)
%
% Arguments:
%   aaa
%
% Contact:
%   jules.van.irsel.gr@dartmouth.edu
%
% Revisions:
%   07/23/2024  initial implementation (jvi)
%

function [bc, resnorm, E2_bg, E3_bg, d2_int, d3_int, weight0, bound] = replicate(track, image, xg, opts)
arguments
    track (1, 1) struct {mustBeNonempty}
    image (1, 1) struct {mustBeNonempty}
    xg (1, 1) struct {mustBeNonempty}
    opts.driver {mustBeMember(opts.driver, ["current", "flow"])} = "flow"
    opts.boundaries (2, 2, :) double = nan(2, 2, 1)
    opts.upsample (1, 1) int16 {mustBePositive} = 1
    opts.num_replications (1, 1) int32 {mustBePositive} = 256
    opts.pos_type {mustBeMember(opts.pos_type, ["angular", "linear"])} = "angular"
    opts.flow_bg (1, 2) double {mustBeNumeric} = nan(1, 2)
    opts.data_smoothing_window (1, 1) int32 {mustBeNumeric} = -1
    opts.do_rotate (1, 1) logical = true
    opts.do_scale (1, 1) logical = true
    opts.arc_definition {mustBeMember(opts.arc_definition, ["conductance", "Pedersen", "Hall", "flux"])} = "Pedersen"
    opts.edge_method {mustBeMember(opts.edge_method, ["contour", "sobel"])} = "contour"
    opts.contour_values (1, 2) double = nan(1, 2)
    opts.boundary_smoothing_window (1, 1) int32 {mustBePositive} = 1
    opts.swap_primary (1, 1) logical = false
    opts.weighting_scale_length (1, 1) double {mustBePositive} = 2e5
    opts.fit_harmonic (1, 1) logical = true
    opts.harmonic_mask (1, 3) double {mustBeNonnegative} = [1, 1, 2] * 5e3
    opts.add_phi_background (1, 1) logical = false
    opts.plot_bg (1, 1) logical = false
    opts.show_plots (1, 1) logical = false
    opts.save_plots (1, 4) logical = false(1, 5)
    opts.save_data (1, 1) logical = false
    opts.auto_lim (1, 1) logical = true
    opts.direc_out (1, :) char {mustBeFolder} = 'output'
    opts.suffix (1, :) char = ''
    opts.starting_letter (1, 1) char = 'A'
end

%% initialize
ftn = 'Arial';
fts = 10 * 2;
lw = 1.4;
qnt = 0.97;
scs = get(0, 'ScreenSize');
% tpl = [0.0 * scs(3), 0.525 * scs(4), 0.5 * scs(3), 0.4 * scs(4)];
% tpr = [0.5 * scs(3), 0.525 * scs(4), 0.5 * scs(3), 0.4 * scs(4)];
% btl = [0.0 * scs(3), 0.050 * scs(4), 0.5 * scs(3), 0.4 * scs(4)];
% btr = [0.5 * scs(3), 0.050 * scs(4), 0.5 * scs(3), 0.4 * scs(4)];
ful = [0.0 * scs(3), 0.050 * scs(4), 1.0 * scs(3), 0.875 * scs(4)];

close all
reset(0)
set(0, 'defaultFigurePaperUnits', 'inches')
set(0, 'defaultTiledlayoutPadding', 'tight')
set(0, 'defaultTiledlayoutTileSpacing', 'tight')
set(0, 'defaultSurfaceEdgeColor', 'flat')
set(0, 'defaultLineLineWidth', lw)
set(0, 'defaultScatterLineWidth', lw)
set(0, 'defaultQuiverLineWidth', lw * 0.7)
aurogem.tools.setall(0, 'FontName', ftn)
aurogem.tools.setall(0, 'FontSize', fts)
aurogem.tools.setall(0, 'Multiplier', 1)

colorcet = @aurogem.tools.colorcet;

scl.x = 1e-3; scl.v = 1e-3; scl.dv = 1e3; scl.p = 1e-3; scl.j = 1e6;
unt.x = 'km'; unt.v = 'km/s'; unt.dv = 'mHz'; unt.p = 'kV'; unt.j = 'uA/m^2';
clm.v = 'D2'; clm.dv = 'CBD1'; clm.p = 'D10'; clm.j = 'D1A'; clm.w = 'L6';
lim.x = [-1, 1] * 125; lim.y = [-1, 1] * 59;  lim.v = [-1, 1] * 1.3; lim.dv = [-1, 1] * 0.1;

scl.vec = 12 * scl.v;

lbl.x = sprintf('Mag. E (%s)', unt.x);
lbl.y = sprintf('Mag. N (%s)', unt.x);
lbl.vx = sprintf('v_E (%s)', unt.v);
lbl.vy = sprintf('v_N (%s)', unt.v);
lbl.j = sprintf('j_{||} (%s)', unt.j);

fid = fopen(fullfile(opts.direc_out, 'output.txt'), 'w');
if not(isempty(opts.suffix))
    opts.suffix = ['_', opts.suffix];
end
if strcmp(opts.driver, 'flow')
    flow_driven = true;
    scl.d = scl.v;
    unt.d = unt.v;
elseif strcmp(opts.driver, 'current')
    flow_driven = false;
    opts.add_phi_background = false;
    opts.do_rotate = false;
    opts.fit_harmonic = false;
    opts.plot_bg = false;
    scl.d = scl.j;
    unt.d = unt.j;
end

%% unpack grid
if opts.upsample > 1
    n_us = opts.upsample;
    for i=[1, fid]; fprintf(i, 'Upsampling grid by a factor of %i.\n', n_us); end
    theta_tmp = linspace(max(xg.theta(1, 1, :)), min(xg.theta(1, 1, :)), xg.lx(3) * n_us);
    phi_tmp = linspace(min(xg.phi(1, :, 1)), max(xg.phi(1, :, 1)), xg.lx(2) * n_us);
    xg_us.theta = permute(repmat(theta_tmp, [xg.lx(1), 1, xg.lx(2) * n_us]), [1, 3, 2]);
    xg_us.phi = permute(repmat(phi_tmp, [xg.lx(1), 1, xg.lx(3) * n_us]), [1, 2, 3]);
    xg_us.x2 = linspace(min(xg.x2), max(xg.x2), (length(xg.x2) - 4) * n_us + 4);
    xg_us.x3 = linspace(min(xg.x3), max(xg.x3), (length(xg.x3) - 4) * n_us + 4);
    xg_us.dx2h = ones(1, numel(xg.dx2h) * n_us) * mean(xg.dx2h);
    xg_us.dx3h = ones(1, numel(xg.dx3h) * n_us) * mean(xg.dx3h);
    xg_us.lx = [xg.lx(1), xg.lx(2) * n_us, xg.lx(3) * n_us];
    xg_us.Bmag = xg.Bmag;
    xg = xg_us;
end

mlat = 90 - squeeze(xg.theta(1, 1, :)) * 180 / pi;
mlon = squeeze(xg.phi(1, :, 1)) * 180 / pi;
x2 = double(xg.x2(3:end - 2));
x3 = double(xg.x3(3:end - 2));
if size(x2, 1)~=1; x2=x2'; end
if size(x3, 1)~=1; x3=x3'; end
dx2 = xg.dx2h;
dx3 = xg.dx3h;
lx2 = xg.lx(2); lx3 = xg.lx(3);
[X2, X3] = ndgrid(x2, x3);
mlon_to_x2 = griddedInterpolant(mlon, x2);
mlat_to_x3 = griddedInterpolant(mlat, x3);
Bmag = abs(mean(xg.Bmag, 'all'));
ar = [range(x2), range(x3), range(x3)];

if opts.auto_lim
    lim.x = [-1, 1] * max(x2) * scl.x;
    lim.y = [-1, 1] * max(x3) * scl.x;
end

%% unpack image
if isfield(image, 'pos_type')
    opts.pos_type = image.pos_type;
end
if strcmp(opts.pos_type, "angular")
    mlon_imag = image.pos(:, :, 1);
    mlat_imag = image.pos(:, :, 2);
    X2_imag = mlon_to_x2(mlon_imag);
    X3_imag = mlat_to_x3(mlat_imag);
else
    X2_imag = image.pos(:, :, 1);
    X3_imag = image.pos(:, :, 2);
end

% assertions
assert(isequal(size(image.pos, 1:2), size(image.flux)), ...
    'image.flux and image.pos(:, :, 1) must have equal sizes.')
assert(size(image.pos, 3)==2, ...
    'Third dimension of image.pos must be 2.')
% assert(issorted(opts.contour_values), ...
%     'Contour values must be in ascending order.')

Q = image.flux;
E0 = image.energy;
if all(isfield(image, 'pedersen'), isfield(image, 'hall'))
    for i=[1, fid]; fprintf(i, 'Pedersen and Hall conductances found.\n'); end
    SIGP = image.pedersen;
    SIGH = image.hall;
else
    Ebar = E0 / 1e3;
    SIGP = 40 * Ebar .* sqrt(Q) ./ (16 + Ebar.^2); % Robinson et al. (1987), Eq. (3)
    SIGH = 0.45 * (Ebar^0.85); % Robinson et al. (1987), Eq. (4)
end
x2_imag = X2_imag(:, 1)';
x3_imag = X3_imag(1, :);

if strcmp(opts.arc_definition, 'conductance') || strcmp(opts.arc_definition, 'Pedersen')
    arc = SIGP.^2;
    ap = 1/2;
    scl.arc = 1e0;
    unt.arc = 'S';
    clm.arc = 'L18';
    lbl.arc = sprintf('Pedersen conductance (%s)', unt.arc);
elseif strcmp(opts.arc_definition, 'Hall')
    arc = SIGH.^2;
    ap = 1/2;
    scl.arc = 1e0;
    unt.arc = 'S';
    clm.arc = 'L18';
    lbl.arc = sprintf('Hall conductance (%s)', unt.arc);
elseif strcmp(opts.arc_definition, 'flux')
    arc = Q;
    ap = 1;
    scl.arc = 1e0;
    unt.arc = 'mW/m^2';
    clm.arc = 'L19';
    lbl.arc = sprintf('Total energy flux (%s)', unt.arc);
end

%% calculate boundaries
if all(isnan(opts.boundaries(:)))
    num_bounds = 2;
    edges = aurogem.tools.find_max_edges(arc, theta=0);
    bsw = opts.boundary_smoothing_window;
    for i=[1, fid]; fprintf(i, 'Boundary smoothing window is approximately %.0f meters.\n', ...
            mean(dx3) * double(bsw)); end
    if strcmp(opts.edge_method, 'sobel')
        bounds = nan(num_bounds, size(edges, 1));
        for i = 1:size(edges, 1)
            [~, bounds(:, i)] = aurogem.tools.peak_detect(edges(i, :), ...
                num=num_bounds, smoothness=0.009);
        end
        x2_bounds_A = sort(x2_imag(2:end - 1));
        x3_bounds_A = smoothdata(x3_imag(bounds(1, :) + 1), 2, 'gaussian', bsw)';
        x2_bounds_B = x2_bounds_A;
        x3_bounds_B = smoothdata(x3_imag(bounds(2, :) + 1), 2, 'gaussian', bsw)';
    elseif strcmp(opts.edge_method, 'contour')
        if all(isnan(opts.contour_values))
            edge_id = round(size(edges, 1) / 2);
            [~, cntr_ids] = aurogem.tools.peak_detect(edges(edge_id, :), ...
                num=num_bounds, smoothness=0.009);
            cntr_vals = arc(edge_id, cntr_ids);
        else
            cntr_vals = opts.contour_values.^(1/ap);
        end
        cntr_A = contour(X2_imag, X3_imag, arc, [1, 1] * cntr_vals(1));
        cntr_B = contour(X2_imag, X3_imag, arc, [1, 1] * cntr_vals(2));
        close(gcf)
        for i=[1, fid]; fprintf(i, 'Primary contour line at %.2f %s.\n', ...
                cntr_A(1, 1)^ap * scl.arc, unt.arc); end
        for i=[1, fid]; fprintf(i, 'Secondary contour line at %.2f %s.\n', ...
                cntr_B(1, 1)^ap * scl.arc, unt.arc); end
        % min_sep = range(x2_imag) * 0.9;
        try
            [x2_bounds_A, x3_bounds_A] = aurogem.tools.unpack_contours(cntr_A, rank=1);
            [x2_bounds_B, x3_bounds_B] = aurogem.tools.unpack_contours(cntr_B, rank=-1);
        catch
            figure
            hold on
            contour(X2_imag, X3_imag, arc, 20)
            contour(X2_imag, X3_imag, arc, [1, 1] * cntr_vals(1), '--r')
            filename = fullfile(opts.direc_out, 'contour_dump.png');
            saveas(gcf, filename)
            error('On or more contours not found. See %s for details.', filename)
        end
    
    
        if min(x2_bounds_A) > min(X2_imag(:)) || max(x2_bounds_A) < max(X2_imag(:))
            warning('Primary boundary does not span image space.')
        end
        if min(x2_bounds_B) > min(X2_imag(:)) || max(x2_bounds_B) < max(X2_imag(:))
            warning('Secondary boundary does not span image space.')
        end
    
        x2_bounds_A = smoothdata(x2_bounds_A, "gaussian");
        x3_bounds_A = smoothdata(x3_bounds_A, "gaussian", bsw);
        x2_bounds_B = smoothdata(x2_bounds_B, "gaussian");
        x3_bounds_B = smoothdata(x3_bounds_B, "gaussian", bsw);
        [~, sort_ids_A] = sort(x2_bounds_A);
        [~, sort_ids_B] = sort(x2_bounds_B);
        x2_bounds_A = aurogem.tools.minsmooth(x2_bounds_A(sort_ids_A));
        x3_bounds_A = x3_bounds_A(sort_ids_A);
        x2_bounds_B = aurogem.tools.minsmooth(x2_bounds_B(sort_ids_B));
        x3_bounds_B = x3_bounds_B(sort_ids_B);
    end
    
    bound.A = griddedInterpolant(x2_bounds_A, x3_bounds_A);
    bound.B = griddedInterpolant(x2_bounds_B, x3_bounds_B);
    if opts.swap_primary
        bound_A_tmp = bound.A;
        bound.A = bound.B;
        bound.B = bound_A_tmp;
    end
    if opts.save_data
        save(fullfile('data', 'boundaries.mat'), 'bound')
    end
else
    x2_bounds_A = squeeze(opts.boundaries(1, 1, :));
    x3_bounds_A = squeeze(opts.boundaries(1, 2, :));
    x2_bounds_B = squeeze(opts.boundaries(2, 1, :));
    x3_bounds_B = squeeze(opts.boundaries(2, 2, :));
    bound.A = griddedInterpolant(x2_bounds_A, x3_bounds_A);
    bound.B = griddedInterpolant(x2_bounds_B, x3_bounds_B);
end

bound_pts = linspace(min(x2_imag) * 1.1, max(x2_imag) * 1.1, length(x2_imag));
angle = griddedInterpolant(bound_pts(1:end - 1), ...
    atan2(smoothdata(diff(bound.A(bound_pts)), "loess") ...
    , diff(bound_pts)));

if opts.show_plots || opts.save_plots(1)
    figure
    set(gcf, 'PaperPosition', [0, 0, 13.2, 5], 'Position', ful)
    t = tiledlayout(1, 2);
    t.Title.String = 'Boundary Definitions';
    t.Title.FontSize = 1.5 * fts;
    ltr = opts.starting_letter;

    nexttile
    text(0.04 - 0.01, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2_imag * scl.x, X3_imag * scl.x, arc.^ap * scl.arc)
    plot(bound_pts * scl.x, bound.A(bound_pts) * scl.x, 'k')
    plot(bound_pts * scl.x, bound.B(bound_pts) * scl.x, '--k')
    xlim(lim.x); ylim(lim.y)
    xlabel(lbl.x); ylabel(lbl.y)
    colormap(gca, colorcet(clm.arc))
    pbaspect(ar)
    legend('', 'Primary', 'Secondary' ...
        , 'Location', 'northwest', 'Orientation', 'horizontal')

    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8)
    hold on
    [C, h] = contour(X2_imag * scl.x, X3_imag * scl.x, arc.^ap * scl.arc, 10, 'ShowText', 1, 'LabelFormat', '%0.1f');
    clabel(C, h, 'Fontsize', 5)
    plot(bound_pts * scl.x, bound.A(bound_pts) * scl.x, 'k')
    plot(bound_pts * scl.x, bound.B(bound_pts) * scl.x, '--k')
    xlim(lim.x); ylim(lim.y)
    yticks([])
    xlabel(lbl.x)
    colormap(gca, colorcet(clm.arc))
    clb = colorbar;
    clb.Label.String = lbl.arc;
    pbaspect(ar)

    if opts.save_plots(1)
        filename = sprintf('boundary_definitions%s.png', opts.suffix);
        filename = fullfile(opts.direc_out, filename);
        fprintf('Saving %s.\n', filename)
        saveas(gcf, filename)
    end
end

%% unpack track(s)
if isstruct(struct2array(track))
    tracks = track;
    track_names = fieldnames(tracks);
    num_tracks = length(track_names);
    if num_tracks == 1
        track = tracks.(cell2mat(track_names(1)));
    else
        for i=[1, fid]; fprintf(i, 'Multiple tracks.\n'); end
    end
else
    num_tracks = 1;
end

x2_traj_p = [];
x3_traj_p = [];
v2_traj_p = [];
v3_traj_p = [];
mask_tracks = false(size(X2));

for track_id = 1:num_tracks
    if num_tracks ~=1
        track_name = cell2mat(track_names(track_id));
        track = tracks.(track_name);
        for i=[1, fid]; fprintf(i, [pad(sprintf(' Current track = %s ', track_name) ...
                , 80, 'both', '-'), '\n']); end
        title_pfx = sprintf('Track %s: ', track_name);
    else
        title_pfx = '';
    end

    % assertions
    assert(isequal(size(track.pos), size(track.flow)), ...
        'track.pos and track.flow must have equal sizes.')
    assert(size(track.pos, 2)==2, ...
        'Second dimension of track.pos and track.flow must be 2.')
    assert(size(track.fac, 1)==size(track.pos, 1), ...
        'Second dimension of track.pos and track.flow must be 2.')

    % ensure northbound trajectory
    if not(issorted(track.pos(:, 2)))
        track.pos = flip(track.pos);
        track.flow = flip(track.flow);
    end

    % unpack TRACK data
    if strcmp(opts.pos_type, "angular")
        mlon_traj = smoothdata(track.pos(:, 1), "loess", 16);
        mlat_traj = smoothdata(track.pos(:, 2), "loess", 16);
        x2_traj = mlon_to_x2(mlon_traj);
        x3_traj = mlat_to_x3(mlat_traj);
    else
        x2_traj = smoothdata(track.pos(:, 1), "loess", 16);
        x3_traj = smoothdata(track.pos(:, 2), "loess", 16);
    end

    if flow_driven
        d2_traj = fillmissing(track.flow(:, 1), 'nearest');
        d3_traj = fillmissing(track.flow(:, 2), 'nearest');
    else
        d2_traj = fillmissing(track.fac, 'nearest');
        d3_traj = zeros(size(d2_traj));
    end

    % extrapolate flow data
    % n_ext = 50;
    % x3_traj_old = x3_traj;
    % fv2 = scatteredInterpolant(x2_traj, x3_traj, v2_traj, 'natural', 'nearest');
    % fv3 = scatteredInterpolant(x2_traj, x3_traj, v3_traj, 'natural', 'nearest');
    % x2_traj = [...
    %     x2_traj(1) + (x2_traj(2) - x2_traj(1)) * (-n_ext:-1), ...
    %     x2_traj', ...
    %     x2_traj(end) + (x2_traj(end) - x2_traj(end - 1)) * (1:n_ext) ...
    %     ]';
    % x3_traj = [...
    %     x3_traj(1) + (x3_traj(2) - x3_traj(1)) * (-n_ext:-1), ...
    %     x3_traj', ...
    %     x3_traj(end) + (x3_traj(end) - x3_traj(end - 1)) * (1:n_ext)]';
    % v2_traj = fv2(x2_traj, x3_traj);
    % v3_traj = fv3(x2_traj, x3_traj);

    % adjust x limits
    lim.x(1) = min([x2_traj(x3_traj * scl.x > lim.y(1) & x3_traj * scl.x < lim.y(2)) * scl.x; lim.x(1)]);
    lim.x(2) = max([x2_traj(x3_traj * scl.x > lim.y(1) & x3_traj * scl.x < lim.y(2)) * scl.x; lim.x(2)]);

    %% replicate in situ flow data
    % 0 = original, 1 = replicated, a = primary boundary, b = secondary boundary
    % minmooth needed for unique independent input of griddedinterpolant

    [~, sort_ids] = sort(x2_traj); % used for griddedInterpolant

    % position where original trajectory meets primary boundary
    traj0 = griddedInterpolant(aurogem.tools.minsmooth(x2_traj(sort_ids)), ...
        x3_traj(sort_ids));

    x0a = fzero(@(x)(traj0(x) - bound.A(x)), 0);
    y0a = traj0(x0a);

    % position where original trajectory meets secondary boundary
    x0b = fzero(@(x)(traj0(x) - bound.B(x)), 0);
    y0b = traj0(x0b);
    width0 = sqrt((x0b - x0a)^2 + (y0b - y0a)^2);

    % remove background flow
    if flow_driven
        fd2_traj = griddedInterpolant(x3_traj, d2_traj);
        fd3_traj = griddedInterpolant(x3_traj, d3_traj);
        d2_traj_0a = fd2_traj(y0a);
        d3_traj_0a = fd3_traj(y0a);
        chi0 = angle(x0a);
        if track_id == 1 % only define one background flow
            if all(isnan(opts.flow_bg))
                alpha0 = atan2(d3_traj_0a, d2_traj_0a);
                gamma0 = pi + chi0 - alpha0;
                d2_traj_0a_rot = cos(gamma0) * d2_traj_0a - sin(gamma0) * d3_traj_0a;
                d3_traj_0a_rot = sin(gamma0) * d2_traj_0a + cos(gamma0) * d3_traj_0a;
                v_bg = [d2_traj_0a - d2_traj_0a_rot, d3_traj_0a - d3_traj_0a_rot];
            else
                v_bg = opts.flow_bg;
            end
            for i=[1, fid]
                fprintf(i, 'Background flow %.2f m/s east and %.2f m/s north.\n', v_bg(1), v_bg(2));
            end
        end
        d2_traj = d2_traj - v_bg(1);
        d3_traj = d3_traj - v_bg(2);
    else
        v_bg = [0, 0];
    end

    % smooth flow data
    if opts.data_smoothing_window > 0
        dx_traj = median(sqrt(diff(x2_traj).^2 + diff(x3_traj).^2));
        sample_freq = 1 / dx_traj;
        if track_id == 1
            fsm = opts.data_smoothing_window;
            pass_freq = sample_freq / double(fsm);
        end
        for i=[1, fid]; fprintf(i, 'Data smoothing window is approximately %.0f meters.\n', 1 / pass_freq); end
        d2_traj = smoothdata(d2_traj, 'gaussian', round(sample_freq / pass_freq));
        d3_traj = smoothdata(d3_traj, 'gaussian', round(sample_freq / pass_freq));
    end
    % fv2_traj = griddedInterpolant(x3_traj, v2_traj, 'spline');
    % fv3_traj = griddedInterpolant(x3_traj, v3_traj, 'spline');
    % x2_traj = linspace(min(x2_traj), max(x2_traj), traj_us * numel(x2_traj));
    % x3_traj = linspace(min(x3_traj), max(x3_traj), traj_us * numel(x3_traj));
    % v2_traj = fv2_traj(x3_traj);
    % v3_traj = fv3_traj(x3_traj);
    % [~, sort_ids] = sort(x2_traj);
    
    

    % replicate
    dx_min = (x2(1) - max(x2_traj)) * 1.1;
    dx_max = (x2(end) - min(x2_traj)) * 1.1;
    dxs = linspace(dx_min, dx_max, opts.num_replications); % eastward displacements
    x2_traj_rep = nan(length(dxs), length(x2_traj));
    x3_traj_rep = nan(length(dxs), length(x3_traj));
    d2_traj_rep = nan(length(dxs), length(d2_traj));
    d3_traj_rep = nan(length(dxs), length(d3_traj));
    for i = 1:length(dxs)
        % translate
        dx = dxs(i);
        dy = bound.A(x0a + dx) - y0a;
        x2_traj_tra = x2_traj + dx;
        x3_traj_tra = x3_traj + dy;

        % determine width at position 1
        % position where replicated trajectory meets primary boundary
        traj1 = griddedInterpolant(aurogem.tools.minsmooth(x2_traj_tra(sort_ids)), ...
            x3_traj_tra(sort_ids));
        x1a = fzero(@(x)(traj1(x) - bound.A(x)), 0);
        y1a = traj1(x1a);
        % position where replicated trajectory meets secondary boundary
        x1b = fzero(@(x)(traj1(x) - bound.B(x)), 0);
        y1b = traj1(x1b);
        width1 = sqrt((x1b - x1a)^2 + (y1b - y1a)^2);
        beta1 = atan2(x1b - x1a, y1b - y1a); % angle b/w line1 and vertical

        if opts.do_scale
            % rotate about (x1a, y1a) by beta
            x2_traj_rot = cos(beta1) * (x2_traj_tra - x1a) ...
                - sin(beta1) * (x3_traj_tra - y1a) + x1a;
            x3_traj_rot = sin(beta1) * (x2_traj_tra - x1a) ...
                + cos(beta1) * (x3_traj_tra - y1a) + y1a;

            % scale about p1a
            scale = width1 / width0;
            x2_traj_scl = x2_traj_rot;
            x3_traj_scl = scale * (x3_traj_rot - y1a) + y1a;

            % rotate back
            x2_traj_rep(i, :) = cos(-beta1) * (x2_traj_scl - x1a) ...
                - sin(-beta1) * (x3_traj_scl - y1a) + x1a;
            x3_traj_rep(i, :) = sin(-beta1) * (x2_traj_scl - x1a) ...
                + cos(-beta1) * (x3_traj_scl - y1a) + y1a;
        else
            x2_traj_rep(i, :) = x2_traj_tra;
            x3_traj_rep(i, :) = x3_traj_tra;
        end

        if opts.do_rotate
            % rotate flows to be tangent to bound_a
            chi1 = angle(x1a) - chi0;
            d2_traj_rep(i, :) = cos(chi1) * d2_traj - sin(chi1) * d3_traj;
            d3_traj_rep(i, :) = sin(chi1) * d2_traj + cos(chi1) * d3_traj;
        else
            d2_traj_rep(i, :) = d2_traj;
            d3_traj_rep(i, :) = d3_traj;
        end

        if i == round(opts.num_replications * 0.3)
            i_p1 = i;
            x2_traj_tra_p1 = x2_traj_tra; x3_traj_tra_p1 = x3_traj_tra;
            x1a_p1 = x1a; y1a_p1 = y1a;
        end
        if i == round(opts.num_replications * 0.7)
            i_p2 = i;
            x2_traj_tra_p2 = x2_traj_tra; x3_traj_tra_p2 = x3_traj_tra;
            x1a_p2 = x1a; y1a_p2 = y1a;
        end
    end

    if opts.save_data
        save(fullfile('data', 'reps.mat'), 'x2_traj_rep', 'x3_traj_rep', 'd2_traj_rep', 'd3_traj_rep')
    end

    if opts.show_plots || opts.save_plots(2)
        figure
        set(gcf, 'PaperPosition', [0, 0, 13.2, 5], 'Position', ful)
        t = tiledlayout(1, 2);
        t.Title.String = sprintf('%sReplication Subset', title_pfx);
        t.Title.FontSize = 1.5 * fts;
        ltr = opts.starting_letter;

        [~, j_p] = min(abs(x3_traj - bound.B(x0b)));

        nexttile
        ms = 400;
        text(0.04 - 0.01, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
        hold on
        pcolor(X2_imag * scl.x, X3_imag * scl.x, arc.^ap * scl.arc)
        plot(bound_pts * scl.x, bound.A(bound_pts) * scl.x, 'k')
        plot(bound_pts * scl.x, bound.B(bound_pts) * scl.x, '--k')
        plot(x2_traj * scl.x, x3_traj * scl.x, 'r')
        plot(x2_traj_tra_p1 * scl.x, x3_traj_tra_p1 * scl.x, '--r')
        plot(x2_traj_tra_p2 * scl.x, x3_traj_tra_p2 * scl.x, '--r')
        scatter(x2_traj_tra_p1(j_p) * scl.x, x3_traj_tra_p1(j_p) * scl.x, ms, 'xr')
        scatter(x2_traj_tra_p2(j_p) * scl.x, x3_traj_tra_p2(j_p) * scl.x, ms, 'xr')
        plot(x2_traj_rep(i_p1, :) * scl.x, x3_traj_rep(i_p1, :) * scl.x, 'b')
        plot(x2_traj_rep(i_p2, :) * scl.x, x3_traj_rep(i_p2, :) * scl.x, 'b')
        scatter(x2_traj_rep(i_p1, j_p) * scl.x, x3_traj_rep(i_p1, j_p) * scl.x, ms, 'xb')
        scatter(x2_traj_rep(i_p2, j_p) * scl.x, x3_traj_rep(i_p2, j_p) * scl.x, ms, 'xb')
        scatter(x0a * scl.x, y0a * scl.x, ms, 'xk')
        scatter(x0b * scl.x, y0b * scl.x, ms, 'xr')
        scatter(x1a_p1 * scl.x, y1a_p1 * scl.x, ms, 'xk')
        scatter(x1a_p2 * scl.x, y1a_p2 * scl.x, ms, 'xk')
        xlim(lim.x); ylim(lim.y)
        xlabel(lbl.x); ylabel(lbl.y)
        colormap(gca, colorcet(clm.arc))
        pbaspect(ar)

        nexttile
        vs = 10;
        cad = 16;
        text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8)
        hold on
        pcolor(X2_imag * scl.x, X3_imag * scl.x, arc.^ap * scl.arc)
        plot(bound_pts * scl.x, bound.A(bound_pts) * scl.x, 'k')
        plot(bound_pts * scl.x, bound.B(bound_pts) * scl.x, '--k')
        quiver(x2_traj_rep(1:cad:end) * scl.x, x3_traj_rep(1:cad:end) * scl.x, ...
            d2_traj_rep(1:cad:end) * scl.d * vs, d3_traj_rep(1:cad:end) * scl.d * vs, 0, '.-b')
        quiver(x2_traj * scl.x, x3_traj * scl.x, ...
            d2_traj * scl.d * vs, d3_traj * scl.d * vs, 0, '.-r')
        xlim(lim.x); ylim(lim.y)
        yticks([])
        xlabel(lbl.x)
        colormap(gca, colorcet(clm.arc))
        pbaspect(ar)
        rectangle('Position', [lim.x(1) + range(lim.x) * 0.7, lim.y(1), range(lim.x) * 0.3, ...
            range(lim.y) / 8], 'FaceColor', 'w', 'FaceAlpha', 0.7, 'EdgeColor', 'none')
        quiver(lim.x(1) + 1.03 * range(lim.x) * 0.7, lim.y(1) + range(lim.y) / 16, ...
            vs, 0, 0, '.-r', 'LineWidth', lw)
        text(lim.x(1) + 1.12 * range(lim.x) * 0.7, lim.y(1) + range(lim.y) / 16, ...
            sprintf('%.0f %s', 1, unt.d), 'FontSize', fts)

        if opts.save_plots(2)
            if num_tracks == 1
                filename = sprintf('replication_subset%s.png', opts.suffix);
            else
                filename = sprintf('replication_subset_%s%s.png' ...
                    , track_name, opts.suffix);
            end
            filename = fullfile(opts.direc_out, filename);
            fprintf('Saving %s.\n', filename)
            saveas(gcf, filename)
        end
    end

    %% interpolate replicated flow data
    d2_int = griddata(x2_traj_rep(:), x3_traj_rep(:), d2_traj_rep(:), X2, X3, 'cubic'); %#ok<*GRIDD> need up to cubic
    d3_int = griddata(x2_traj_rep(:), x3_traj_rep(:), d3_traj_rep(:), X2, X3, 'cubic');
    fprintf('Done interpolating replicated data.\n')

    if num_tracks > 1
        interpolations.(track_name).d2_int = d2_int;
        interpolations.(track_name).d3_int = d3_int;
    end
    x2_traj_p = [x2_traj_p; x2_traj]; %#ok<AGROW>
    x3_traj_p = [x3_traj_p; x3_traj]; %#ok<AGROW>
    v2_traj_p = [v2_traj_p; d2_traj]; %#ok<AGROW>
    v3_traj_p = [v3_traj_p; d3_traj]; %#ok<AGROW>

    mask_track = false(size(X2)); % select data around trajectory
    for i = 1:lx3
        b = fzero(@(x) traj0(x) - x3(i), 0);
        db = opts.harmonic_mask(3);
        mask_track(:, i) = (x2' >= b - db & x2' < b + db);
    end
    mask_tracks = mask_tracks | mask_track;
end

%% weight function
if num_tracks > 1
    assert(num_tracks==2, '3 or more tracks TBD.\n')

    track_name0 = cell2mat(track_names(1));
    track_name1 = cell2mat(track_names(2));
    x20 = mlon_to_x2(tracks.(track_name0).pos(:, 1));
    x30 = mlat_to_x3(tracks.(track_name0).pos(:, 2));
    x21 = mlon_to_x2(tracks.(track_name1).pos(:, 1));
    x31 = mlat_to_x3(tracks.(track_name1).pos(:, 2));
    l0 = length(x20);
    l1 = length(x21);

    distance0 = nan([size(X2), l0]);
    for i=1:l0
        distance0(:, :, i) = sqrt((X2 - x20(i)).^2 + (X3 - x30(i)).^2);
    end
    distance1 = nan([size(X2), l1]);
    for i=1:l1
        distance1(:, :, i) = sqrt((X2 - x21(i)).^2 + (X3 - x31(i)).^2);
    end

    min_dist0 = min(distance0, [], 3);
    min_dist1 = min(distance1, [], 3);

    weight_scale = opts.weighting_scale_length;
    weight0 = (1 + tanh((min_dist1 - min_dist0) * 2 / weight_scale)) / 2;
    weight0 = aurogem.tools.smoothdata2(weight0, "gaussian", 100);
    weight1 = 1 - weight0;
    for i=[1, fid]; fprintf(i, 'Maximum weight = %f.\n', max(weight0(:))); end

    d2_int0 = interpolations.(track_name0).d2_int;
    d3_int0 = interpolations.(track_name0).d3_int;
    d2_int1 = interpolations.(track_name1).d2_int;
    d3_int1 = interpolations.(track_name1).d3_int;
    d2_int = d2_int0 .* weight0 + d2_int1 .* weight1;
    d3_int = d3_int0 .* weight0 + d3_int1 .* weight1;
else
    weight0 = ones(size(X2));
end

if (opts.show_plots || opts.save_plots(3)) && not(flow_driven)
    max_j = quantile(abs(d2_int(:)), qnt);
    lim.j = [-1, 1] * max_j * scl.j;

    figure
    set(gcf, 'PaperPosition', [0, 0, 13.2, 6.4], 'Position', ful)
    t = tiledlayout(1, 2);
    t.Title.String = sprintf('%sReplication Summary', title_pfx);
    t.Title.FontSize = 1.5 * fts;
    ltr = opts.starting_letter;

    nexttile
    vs = 10;
    text(0.04 - 0.01, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2_imag * scl.x, X3_imag * scl.x, arc.^ap * scl.arc)
    plot(bound_pts * scl.x, bound.A(bound_pts) * scl.x, 'k')
    plot(bound_pts * scl.x, bound.B(bound_pts) * scl.x, '--k')
    quiver(x2_traj * scl.x, x3_traj * scl.x, ...
        d2_traj * scl.d * vs, d3_traj * scl.d * vs, 0, '.-r')
    xlim(lim.x); ylim(lim.y)
    xlabel(lbl.x); ylabel(lbl.y)
    colormap(gca, colorcet(clm.arc))
    clb = colorbar;
    clb.Label.String = lbl.arc;
    clb.Location = 'southoutside';
    pbaspect(ar)

    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8)
    hold on
    pcolor(X2_imag * scl.x, X3_imag * scl.x, d2_int * scl.j)
    plot(bound_pts * scl.x, bound.A(bound_pts) * scl.x, 'k')
    plot(bound_pts * scl.x, bound.B(bound_pts) * scl.x, '--k')
    quiver(x2_traj * scl.x, x3_traj * scl.x, ...
        d2_traj * scl.d * vs, d3_traj * scl.d * vs, 0, '.-r')
    xlim(lim.x); ylim(lim.y)
    yticks([])
    xlabel(lbl.x)
    colormap(gca, colorcet(clm.j))
    clb = colorbar;
    clb.Label.String = lbl.j;
    clb.Location = 'southoutside';
    clim(lim.j)
    pbaspect(ar)
    rectangle('Position', [lim.x(1) + range(lim.x) * 0.7, lim.y(1), range(lim.x) * 0.3, ...
        range(lim.y) / 8], 'FaceColor', 'w', 'FaceAlpha', 0.7, 'EdgeColor', 'none')
    quiver(lim.x(1) + 1.03 * range(lim.x) * 0.7, lim.y(1) + range(lim.y) / 16, ...
        vs, 0, 0, '.-r', 'LineWidth', lw)
    text(lim.x(1) + 1.12 * range(lim.x) * 0.7, lim.y(1) + range(lim.y) / 16, ...
        sprintf('%.0f %s', 1, unt.d), 'FontSize', fts)

    if opts.save_plots(3)
        filename = sprintf('replication_summary%s.png', opts.suffix);
        filename = fullfile(opts.direc_out, filename);
        fprintf('Saving %s.\n', filename)
        saveas(gcf, filename)
    end
end

% if current driven, stop here
if flow_driven
    E2_int =  d3_int * Bmag; % v = ExB/B^2
    E3_int = -d2_int * Bmag;
    E2_bg =  v_bg(2) * Bmag;
    E3_bg = -v_bg(1) * Bmag;
else
    bc = -d2_int; % top boundary condition for gemini is anti-parallel (upward)
    E2_bg =  opts.flow_bg(2) * Bmag;
    E3_bg = -opts.flow_bg(1) * Bmag;
    resnorm = nan(1, 2);
    return
end

%% generate potential map
% translated from Alex Mule's python code Oct 9, 2023
% extrapolate data to avoid edge effects with fast-fourier transforms
nfill = 32;
x2_ext = [x2(1) + dx2(1) * (-nfill:-1), x2, x2(end) + dx2(end) * (1:nfill)];
x3_ext = [x3(1) + dx3(1) * (-nfill:-1), x3, x3(end) + dx3(end) * (1:nfill)];
[X2_ext, X3_ext] = ndgrid(x2_ext, x3_ext);
fE2_int = griddedInterpolant(X2, X3, E2_int, 'spline', 'nearest');
fE3_int = griddedInterpolant(X2, X3, E3_int, 'spline', 'nearest');
E2_int_ext = fE2_int(X2_ext, X3_ext);
E3_int_ext = fE3_int(X2_ext, X3_ext);

% wave vector convention: 2 pi ( -f_Ny : +f_Ny )
k2_Ny = 2 * pi / (2 * mean(dx2));
k3_Ny = 2 * pi / (2 * mean(dx3));
k2 = linspace(-k2_Ny, k2_Ny, lx2 + 2 * nfill);
k3 = linspace(-k3_Ny, k3_Ny, lx3 + 2 * nfill);
[K2, K3] = ndgrid(k2, k3);

% Fourier transform of electric field
G2 = fftshift(fft2(E2_int_ext));
G3 = fftshift(fft2(E3_int_ext));

% Fourier transform of phi
Gphi = 1i * (K2 .* G2 + K3 .* G3) ./ (K2.^2 + K3.^2);

% inverse Fourier transform of Gphi + resampling onto working grid
phi0 = real(ifft2(ifftshift(Gphi)));
phi0 = phi0(nfill + 1:end - nfill, nfill + 1:end - nfill);

% find best fit harmonic function
mask_A = false(size(X2)); % select data around boundary A
mask_B = false(size(X2)); % select data around boundary B
for i = 1:lx2
    b = bound.A(x2(i));
    db = opts.harmonic_mask(1);
    mask_A(i, :) = (x3 >= b - db & x3 < b + db);
end
for i = 1:lx2
    b = bound.B(x2(i));
    db = opts.harmonic_mask(2);
    mask_B(i, :) = (x3 >= b - db & x3 < b + db);
end
mask = mask_A | mask_B | mask_tracks;
xdata = [X2(mask), X3(mask)];
Edata = [E2_int(mask), E3_int(mask)];

if opts.fit_harmonic
    fprintf('Fitting harmonic function.\n')
    [bc, ~] = aurogem.tools.find_harmonic(phi0, xdata, Edata, xg);
else
    [E2_0, E3_0] = gradient(-phi0, mean(dx2), mean(dx3));
    E_0 = mean([E2_0(:);E3_0(:)]) - mean([E2_int(:), E3_int(:)]);
    harm = E_0(1) * X2 + E_0(2) * X3;
    bc = phi0 + harm;
end

% adding back background electric field
if opts.add_phi_background
    bc = bc - X2 * E2_bg - X3 * E3_bg;
end

% remove Nyquist ringing
try
    bc = smoothdata2(bc, "gaussian", 4);
catch
    warning('Potential not smoothed. Need version R2023b or higher. Current version = %s\n', version)
end

% ensure first and second differences match at north and south edges
bc(:, 1) = 3 * bc(:, 3) - 2 * bc(:, 4);
bc(:, 2) = 2 * bc(:, 3) - 1 * bc(:, 4);
bc(:, end - 1) = 2 * bc(:, end - 2) - 1 * bc(:, end - 3);
bc(:, end - 0) = 3 * bc(:, end - 2) - 2 * bc(:, end - 3);

% zero average potential
bc = bc - mean(bc, 'all');

%% plot final flow fields
fSIGP = griddedInterpolant(X2_imag, X3_imag, SIGP);
if range(dx2) < 1e-3
    [E2, E3] = gradient(-bc', mean(dx2), mean(dx3));
else
    [E2, E3] = gradient(-bc', dx2, dx3);
end
v2 = -E3' / Bmag;
v3 =  E2' / Bmag;
divv = divergence(X3, X2, v3, v2);
divv_int = divergence(X3, X2, d3_int, d2_int);
divE = divergence(X3, X2, E3', E2');
fac_divE = fSIGP(X2, X3) .* divE;
v2_err = v2 - d2_int;
v3_err = v3 - d3_int;

dist2 = v2_err(mask);
dist3 = v3_err(mask);
v2_err_min = round(mean(dist2)) - round(2 * std(dist2));
v2_err_max = round(mean(dist2)) + round(2 * std(dist2));
v3_err_min = round(mean(dist3)) - round(2 * std(dist3));
v3_err_max = round(mean(dist3)) + round(2 * std(dist3));
resnorm = [ ...
    norm(abs(v2_err(mask))) / norm(abs(d2_int(mask))) ...
    , norm(abs(v3_err(mask))) / norm(abs(d3_int(mask))) ...
    ];
fprintf([pad('Electrostatic enforcement results:', 80, 'both', '-'), '\n'])
for i=[1, fid]; fprintf(i, 'Eastward error range = %i to %i m/s.\n', v2_err_min, v2_err_max); end
for i=[1, fid]; fprintf(i, 'Northward error range = %i to %i m/s.\n', v3_err_min, v3_err_max); end
for i=[1, fid]; fprintf(i, 'Eastward relative norm of the residuals = %.2f.\n', resnorm(1)); end
for i=[1, fid]; fprintf(i, 'Northward relative norm of the residuals = %.2f.\n', resnorm(2)); end
fprintf([pad('', 80, 'both', '-'), '\n'])

v2_int_p = (d2_int + opts.plot_bg * v_bg(1)) * scl.v;
v3_int_p = (d3_int + opts.plot_bg * v_bg(2)) * scl.v;
v2_traj_p = (v2_traj_p + opts.plot_bg * v_bg(1)) * scl.v;
v3_traj_p = (v3_traj_p + opts.plot_bg * v_bg(2)) * scl.v;
v2_traj_pv = v2_traj_p * scl.vec / scl.v;
v3_traj_pv = v3_traj_p * scl.vec / scl.v;
v2_p = (v2 + opts.plot_bg * v_bg(1)) * scl.v;
v3_p = (v3 + opts.plot_bg * v_bg(2)) * scl.v;

if opts.show_plots || opts.save_plots(3)
    figure
    set(gcf, 'PaperPosition', [0, 0, 13.2, 6.4], 'Position', ful)
    tiledlayout(3, 3);
    ltr = opts.starting_letter;

    if opts.auto_lim
        max_v = quantile(abs([v2(:) + v_bg(1);v3(:) + v_bg(2)]), qnt);
        max_dv = quantile(abs([divv(:);divv_int(:)]), qnt);
        max_p = quantile(abs(bc(:)), qnt);
        max_j = quantile(abs(fac_divE(:)), qnt);
        lim.v = [-1, 1] * max_v * scl.v;
        lim.dv = [-1, 1] * max_dv * scl.dv;
        lim.p = [-1, 1] * max_p * scl.p;
        lim.j = [-1, 1] * max_j * scl.j;
    end

    % row 1
    nexttile
    title('Interpolated flow')
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, v2_int_p)
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    colormap(gca, colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = lbl.vx;
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([])
    ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    title('Helmholtz decomposition')
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, v2_p)
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    colormap(gca, colorcet(clm.v))
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([]); yticks([])
    pbaspect(ar)

    nexttile
    title('Difference & Potential')
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, v2_err * scl.v)
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    contour(X2 * scl.x, X3 * scl.x, double(mask), [0.5, 0.5], 'k')
    colormap(gca, colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('\\Delta %s', lbl.vx);
    xlim(lim.x); ylim(lim.y); clim(lim.v / 3)
    xticks([]); yticks([])
    pbaspect(ar)

    %row 2
    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, v3_int_p)
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    colormap(gca, colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = lbl.vy;
    clb.Location = 'westoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([])
    ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, v3_p)
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    colormap(gca, colorcet(clm.v))
    xlim(lim.x); ylim(lim.y); clim(lim.v)
    xticks([]); yticks([])
    pbaspect(ar)

    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, v3_err * scl.v)
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    contour(X2 * scl.x, X3 * scl.x, double(mask), [0.5, 0.5], 'k')
    colormap(gca, colorcet(clm.v))
    clb = colorbar;
    clb.Label.String = sprintf('\\Delta %s', lbl.vy);
    xlim(lim.x); ylim(lim.y); clim(lim.v / 3)
    xticks([]); yticks([])
    pbaspect(ar)

    % row 3
    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, divv_int * scl.dv)
    colormap(gca, colorcet(clm.dv))
    clb = colorbar;
    clb.Label.String = sprintf('\\nabla\\cdot\\bfv\\rm (%s)', unt.dv);
    clb.Location = 'southoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.dv)
    xlabel(lbl.x); ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, fac_divE * scl.j)
    colormap(gca, colorcet(clm.j))
    clb = colorbar;
    clb.Label.String = sprintf('\\Sigma_P\\nabla\\cdot\\bfE\\rm (%s)', unt.j);
    clb.Location = 'southoutside';
    xlim(lim.x); ylim(lim.y); clim(lim.j)
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    nexttile
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8)
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, bc * scl.p)
    colormap(gca, colorcet(clm.p))
    clb = colorbar;
    clb.Label.String = sprintf('Potential (%s)', unt.p);
    clb.Location = 'southoutside';
    xlim(lim.x); ylim(lim.y);
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    if opts.save_plots(3)
        filename = sprintf('replication_summary%s.png', opts.suffix);
        filename = fullfile(opts.direc_out, filename);
        fprintf('Saving %s.\n', filename)
        saveas(gcf, filename)
    end
end

if opts.show_plots || opts.save_plots(4)
    figure
    set(gcf, 'PaperPosition', [0, 0, 13.2, 4.5], 'Position', ful)
    tiledlayout(1, 3);
    ltr = opts.starting_letter;

    % hsv plot variables
    MLAT = 90 - xg.theta * 180 / pi;
    MLON = xg.phi * 180 / pi;
    ALT = xg.alt;
    V2 = permute(repmat(v2 + v_bg(1) * opts.plot_bg, [1, 1, size(MLAT, 1)]), [3, 1, 2]);
    V3 = permute(repmat(v3 + v_bg(2) * opts.plot_bg, [1, 1, size(MLAT, 1)]), [3, 1, 2]);
    V2_err = permute(repmat(v2_err, [1, 1, size(MLAT, 1)]), [3, 1, 2]);
    V3_err = permute(repmat(v3_err, [1, 1, size(MLAT, 1)]), [3, 1, 2]);
    mlon_ref = mean(MLON(:));
    hsv_sat = 1e3;
    eps = 0.02;
    [hsv_map_clb, ~, ~, hsv_alt, hsv_alt_map] = ...
        aurogem.tools.hsv_params(V2, V3, MLAT, MLON, ALT, 300e3, mlon_ref, hsv_sat);
    [hsv_map_clb_err, ~, ~, hsv_alt_err, hsv_alt_map_err] = ...
        aurogem.tools.hsv_params(V2_err, V3_err, MLAT, MLON, ALT, 300e3, mlon_ref, hsv_sat * 0.3);

    if opts.auto_lim
        qnt = 0.99;
        max_v = quantile(abs([v2(:) + v_bg(1);v3(:) + v_bg(2)]), qnt);
        max_dv = quantile(abs([divv(:);divv_int(:)]), qnt);
        max_p = quantile(abs(bc(:)), qnt);
        max_j = quantile(abs(fac_divE(:)), qnt);
        lim.v = [-1, 1] * max_v * scl.v;
        lim.dv = [-1, 1] * max_dv * scl.dv;
        lim.p = [-1, 1] * max_p * scl.p;
        lim.j = [-1, 1] * max_j * scl.j;
    end

    nexttile
    title('Input flow')
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, hsv_alt);
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    quiver(60, -51, 1e3 * scl.vec, 0, 0, '.-r')
    text(75, -50, '1 km/s', 'FontSize', 0.8 * fts, 'Color', 'w')
    contour(X2_imag * scl.x, X3_imag * scl.x, arc.^ap * scl.arc, 4, '--k')
    colormap(gca, hsv_alt_map)
    clb = colorbar;
    colormap(clb, hsv_map_clb)
    clb.Limits = [0, 1];
    clb.Ticks = [0 + eps, 1/4, 1/2, 3/4, 1 - eps];
    clb.TickLabels = {'W', 'S', 'E', 'N', 'W'};
    clb.Label.String = ['Sat. at ', num2str(hsv_sat * scl.v), ' ', unt.v];
    clb.Location = 'southoutside';
    clim([0, 1])
    xlim(lim.x); ylim(lim.y);
    xlabel(lbl.x); ylabel(lbl.y)
    pbaspect(ar)

    nexttile
    title('Input - interpolated flow')
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, hsv_alt_err);
    contour(X2 * scl.x, X3 * scl.x, double(mask), [0.5, 0.5], '--k')
    quiver(x2_traj_p * scl.x, x3_traj_p * scl.x, v2_traj_pv, v3_traj_pv, 0, '.-r')
    colormap(gca, hsv_alt_map_err)
    clb = colorbar;
    colormap(clb, hsv_map_clb_err)
    clb.Limits = [0, 1];
    clb.Ticks = [0 + eps, 1/4, 1/2, 3/4, 1 - eps];
    clb.TickLabels = {'W', 'S', 'E', 'N', 'W'};
    clb.Label.String = ['Sat. at ', num2str(hsv_sat * scl.v * 0.3), ' ', unt.v];
    clb.Location = 'southoutside';
    clim([0, 1])
    xlim(lim.x); ylim(lim.y);
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    nexttile
    title('Input potential')
    text(0.04, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8)
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, bc * scl.p)
    colormap(gca, colorcet(clm.p))
    clb = colorbar;
    clb.Label.String = sprintf('Potential (%s)', unt.p);
    clb.Location = 'southoutside';
    xlim(lim.x); ylim(lim.y);
    yticks([])
    xlabel(lbl.x)
    pbaspect(ar)

    if opts.save_plots(4)
        filename = sprintf('replication_flow%s.png', opts.suffix);
        filename = fullfile(opts.direc_out, filename);
        fprintf('Saving %s.\n', filename)
        saveas(gcf, filename)
    end

end

if (opts.show_plots || opts.save_plots(5)) && num_tracks == 2
    figure
    set(gcf, 'PaperPosition', [0, 0, 13.2, 5], 'Position', ful)
    t = tiledlayout(1, 2);
    t.Title.String = 'Weighting Maps';
    t.Title.FontSize = 1.5 * fts;
    ltr = opts.starting_letter;

    nexttile
    text(0.04 - 0.01, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8); ltr = ltr + 1;
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, weight0)
    xlim(lim.x); ylim(lim.y); clim([0, 1])
    xlabel(lbl.x); ylabel(lbl.y)
    clb = colorbar;
    colormap(gca, colorcet(clm.w))
    clb.Label.String = sprintf('Weight %s', track_name0);
    clb.Location = 'southoutside';
    pbaspect(ar)

    nexttile
    text(0.04 - 0.01, 0.9, char(ltr), 'units', 'normalized', 'FontSize', fts * 0.8)
    hold on
    pcolor(X2 * scl.x, X3 * scl.x, weight1)
    xlim(lim.x); ylim(lim.y); clim([0, 1])
    xlabel(lbl.x); ylabel(lbl.y)
    clb = colorbar;
    colormap(gca, colorcet(clm.w))
    clb.Label.String = sprintf('Weight %s', track_name1);
    clb.Location = 'southoutside';
    pbaspect(ar)

    if opts.save_plots(5)
        filename = sprintf('replication_weights%s.png', opts.suffix);
        filename = fullfile(opts.direc_out, filename);
        fprintf('Saving %s.\n', filename)
        saveas(gcf, filename)
    end
end

if ~opts.show_plots
    close all
end
fclose all;
end
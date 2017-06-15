import matplotlib as mpl
import numpy as np
from plot_fits_files import set_up_axes


def alpha_overlay(C_a, a_a, C_b, a_b=None):
    if a_b is None:
        a_b = np.ones(a_a.shape)
    c_a = np.array([a_a.T] * 3).T * C_a
    c_b = np.array([a_b.T] * 3).T * C_b
    c_out = c_a + ((1 - a_a.T) * c_b.T).T
    return c_out


def alpha_maps(maps, colors=None, vmin=0, vmax=15):
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    iter_cycle = iter(mpl.rcParams['axes.prop_cycle'])
    for mdx, m in enumerate(maps):
        if colors is None:
            c = next(iter_cycle)['color']
        else:
            c = colors[mdx]
        base_color = np.array(mpl.colors.to_rgb(c))
        norm_map = norm(m)
        if mdx == 0:
            background_color = np.ones(3)
        background_color = alpha_overlay(base_color, norm_map, background_color)
    return background_color


def make_alpha_bar(color, vmin=-1, vmax=15):
    # vmin of -1 to make lables line up correctly
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    a_a = norm(range(vmin, vmax))
    C_a = np.array(mpl.colors.to_rgb(color))
    new_cm = alpha_overlay(C_a, a_a, np.ones(3))
    return mpl.colors.ListedColormap(new_cm), norm


def plot_alpha_bar(color, grid, ticks=[]):
    bar, norm = make_alpha_bar(color)
    ax_bar = plt.subplot(grid)
    cb = mpl.colorbar.ColorbarBase(ax_bar, cmap=bar, norm=norm, orientation='vertical', ticks=ticks)
    return ax_bar, cb


def plot_masks(gz3d, grid, colors = ['C1', 'C0', 'C3', 'C2']):
    gs = gridspec.GridSpecFromSubplotSpec(1, 5, width_ratios=[0.84, 0.04, 0.04, 0.04, 0.04], wspace=0.015, subplot_spec=grid)
    maps = [gz3d.bar_mask, gz3d.spiral_mask, gz3d.star_mask, gz3d.center_mask]
    all_mask = alpha_maps(maps, colors)
    # plot image
    ax1 = plt.subplot(gs[0], projection=gz3d.wcs)
    set_up_axes(ax1, color_grid='C7', color_tick='black')
    ax1.imshow(all_mask)
    # plot hexagons
    ax1.add_patch(gz3d.get_hexagon())
    ax1.add_patch(gz3d.get_hexagon(correct_hex=True, edgecolor='C4'))
    # plot center and star ellipses
    center_ellip = gz3d.get_center_ellipse_list()
    for e in center_ellip:
        ax1.add_artist(e)
    star_ellip = gz3d.get_star_ellipse_list()
    for e in star_ellip:
        ax1.add_artist(e)
    # make a legend
    bar_patch = mpatches.Patch(color=colors[0], label='Bar')
    spiral_patch = mpatches.Patch(color=colors[1], label='Spiral')
    star_patch = mpatches.Patch(color=colors[2], label='Star')
    center_patch = mpatches.Patch(color=colors[3], label='Center')
    plt.legend(handles=[bar_patch, spiral_patch, star_patch, center_patch], ncol=4, loc='lower center', mode='expand')
    # make colorbars
    ax_bar, cb_bar = plot_alpha_bar(colors[0], gs[1])
    ax_spiral, cb_spiral = plot_alpha_bar(colors[1], gs[2])
    ax_star, cb_star = plot_alpha_bar(colors[2], gs[3])
    ax_center, cb_center = plot_alpha_bar(colors[3], gs[4])
    tick_labels = np.arange(0, 16)
    tick_locs = tick_labels - 0.5
    cb_center.set_ticks(tick_locs)
    cb_center.set_ticklabels(tick_labels)
    cb_center.set_label('Count')


if __name__ == '__main__':
    from gz3d_fits import gz3d_fits
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import matplotlib.gridspec as gridspec
    file_name = '/Volumes/Work/GZ3D/MPL5_fits/1-167242_127_5679242.fits.gz'
    gz3d = gz3d_fits(file_name)
    fig = plt.figure(1)
    gs = gridspec.GridSpec(1, 1)
    plot_masks(gz3d, gs[0])
    plt.show()

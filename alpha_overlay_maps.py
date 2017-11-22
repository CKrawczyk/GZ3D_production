import matplotlib as mpl
import matplotlib.gridspec as gridspec
import numpy as np
from plot_fits_files import set_up_axes, plot_original
import matplotlib.pyplot as plt
from plot_by_r import zero_theta_line


def alpha_overlay(C_a, a_a, C_b, a_b=None):
    if a_b is None:
        a_b = np.ones(a_a.shape)
    c_a = np.array([a_a.T] * 3).T * C_a
    c_b = np.array([a_b.T] * 3).T * C_b
    c_out = c_a + ((1 - a_a.T) * c_b.T).T
    return c_out


def alpha_maps(maps, colors=None, vmin=0, vmax=15):
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
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
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    a_a = norm(range(vmin, vmax))
    C_a = np.array(mpl.colors.to_rgb(color))
    new_cm = alpha_overlay(C_a, a_a, np.ones(3))
    return mpl.colors.ListedColormap(new_cm), norm


def make_alpha_color(count, color, vmin=1, vmax=15):
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    return mpl.colors.to_rgb(color) + (norm(count), )


def plot_alpha_bar(color, grid, ticks=[]):
    bar, norm = make_alpha_bar(color)
    ax_bar = plt.subplot(grid)
    cb = mpl.colorbar.ColorbarBase(ax_bar, cmap=bar, norm=norm, orientation='vertical', ticks=ticks)
    cb.outline.set_linewidth(0)
    return ax_bar, cb


def plot_masks(gz3d, grid, colors=['C1', 'C0', 'C3', 'C2'], sub_grid_ratio=[0.95, 0.05]):
    gs = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=sub_grid_ratio, subplot_spec=grid, wspace=0.01)
    gs_inner = gridspec.GridSpecFromSubplotSpec(1, 4, wspace=0, subplot_spec=gs[1])
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
    for e, count in zip(center_ellip, gz3d.center_clusters['count']):
        e.set_edgecolor(make_alpha_color(count, 'C2'))
        ax1.add_artist(e)
    star_ellip = gz3d.get_star_ellipse_list()
    for e, count in zip(star_ellip, gz3d.star_clusters['count']):
        e.set_edgecolor(make_alpha_color(count, 'C3'))
        ax1.add_artist(e)
    # plot theta=0 line
    x_theta, y_theta = zero_theta_line(gz3d)
    ax1.plot(x_theta, y_theta, 'C5')
    ax1.annotate(r'$\theta$', xy=(0.2, 0.8), xytext=(0.05, 0.9), arrowprops={'arrowstyle': '->', 'connectionstyle': 'angle3'}, xycoords='axes fraction', textcoords='axes fraction')
    # make a legend
    bar_patch = mpl.patches.Patch(color=colors[0], label='Bar')
    spiral_patch = mpl.patches.Patch(color=colors[1], label='Spiral')
    star_patch = mpl.patches.Patch(color=colors[2], label='Star')
    center_patch = mpl.patches.Patch(color=colors[3], label='Center')
    plt.legend(handles=[bar_patch, spiral_patch, star_patch, center_patch], ncol=2, loc='lower center', mode='expand')
    # make colorbars
    ax_bar, cb_bar = plot_alpha_bar(colors[0], gs_inner[0])
    ax_spiral, cb_spiral = plot_alpha_bar(colors[1], gs_inner[1])
    ax_star, cb_star = plot_alpha_bar(colors[2], gs_inner[2])
    ax_center, cb_center = plot_alpha_bar(colors[3], gs_inner[3])
    ax_center.tick_params(axis=u'both', which=u'both', length=0)
    tick_labels = np.arange(0, 16)
    tick_locs = tick_labels - 0.5
    cb_center.set_ticks(tick_locs)
    cb_center.set_ticklabels(tick_labels)
    cb_center.set_label('Count')
    return ax1


def plot_original_with_mask(gz3d, sub_grid_ratio=[0.9, 0.1], fdx=1):
    fig_width = 12
    fig_height = 4.5
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.05, right=0.94, bottom=0.05, top=0.94, wspace=0.05, hspace=0.3)
    fig = plt.figure(fdx, figsize=(fig_width, fig_height))
    ax_00 = plot_original(gz3d, gs[0, 0], sub_grid_ratio=[0.9, 0.1])
    ax_00.set_title('{0}'.format(gz3d.metadata['MANGAID'][0]))
    ax_10 = plot_masks(gz3d, gs[0, 1], sub_grid_ratio=sub_grid_ratio)
    return fig


if __name__ == '__main__':
    from gz3d_fits import gz3d_fits
    import mpl_style
    import os
    plt.style.use(mpl_style.style1)
    '''
    file_name = '/Volumes/Work/GZ3D/MPL5_fits/1-167242_127_5679242.fits.gz'
    gz3d = gz3d_fits(file_name)
    gz3d.get_bpt()
    # fig = plt.figure(1)
    # gs = gridspec.GridSpec(1, 1)
    # plot_masks(gz3d, gs[0], sub_grid_ratio=[0.9, 0.1])
    fig = plot_original_with_mask(gz3d)
    plt.show()
    '''
    filepath = '/Volumes/Work/GZ3D/MPL5_fits'
    output_folder = '/Users/coleman/Desktop/plots_for_talk/masks'
    id_list = [
        '1-163516_127_5679061',
        '1-135044_91_5682572',
        '1-135078_127_5679767',
        '1-135468_127_5679686',
        '1-210923_127_5679193',
        '1-216958_37_5680828',
        '1-246549_127_5679436',
        '1-37211_127_5679377',
        '1-574355_127_5679349',
        '1-167242_127_5679242'
    ]
    for fdx, mid in enumerate(id_list):
        output_name = os.path.join(output_folder, mid)
        input_name = os.path.join(filepath, mid) + '.fits.gz'
        gz3d = gz3d_fits(input_name)
        gz3d.get_bpt()
        fig = plot_original_with_mask(gz3d, fdx=fdx)
        fig.savefig('{0}_masks.png'.format(output_name))
        plt.close(fig)

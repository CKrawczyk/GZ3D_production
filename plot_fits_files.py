import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_style
import os
import progressbar as pb
from gz3d_fits import gz3d_fits
from stack_spec import get_stacked_spectra

plt.style.use(mpl_style.style1)


def set_up_axes(ax, color_grid='white', color_tick='white'):
    # extract the coordinates from the axes object
    ra = ax.coords['ra']
    dec = ax.coords['dec']
    # add axis labels
    ra.set_axislabel('RA')
    dec.set_axislabel('Dec')
    # rotate tick labels on dec
    dec.ticklabels.set_rotation(90)
    # add a coordinate grid to the image
    ax.coords.grid(color=color_grid, alpha=0.5, linestyle='solid', lw=1.5)
    for coord in [ra, dec]:
        # set the tick formats
        coord.set_major_formatter('d.dd')
        coord.set_ticks(color=color_tick)
        coord.display_minor_ticks(True)


def make_ax(grid, projection, **kwargs):
    ax = plt.subplot(grid, projection=projection)
    set_up_axes(ax, **kwargs)
    return ax


def plot_mask(data, kind, grid, cmap=plt.cm.Purples, sub_grid_ratio=[0.95, 0.05], spaxel_grid=False):
    gs_inner = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=sub_grid_ratio, subplot_spec=grid)
    ax = make_ax(gs_inner[0], data.wcs, color_grid='C0', color_tick='black')
    ax.add_patch(data.get_hexagon())
    ax.add_patch(data.get_hexagon(correct_hex=True, edgecolor='C4'))
    if spaxel_grid:
        v_line, h_line = data.get_spaxel_grid()
        ax.plot(v_line[0], v_line[1], color='C7', lw=0.5, alpha=0.3)
        ax.plot(h_line[0], h_line[1], color='C7', lw=0.5, alpha=0.3)
    image = getattr(data, '{0}_mask'.format(kind))
    vmax = image.max()
    if vmax != 0:
        im = ax.imshow(image, cmap=cmap, vmin=0, vmax=vmax)
        ax_bar = plt.subplot(gs_inner[1])
        cb = plt.colorbar(im, cax=ax_bar)
        cb.set_label('Count')
    else:
        im = ax.imshow(image, cmap=cmap)
    return ax


def plot_ellipse(data, kind, *mask_params, **mask_kwargs):
    ax = plot_mask(data, kind, *mask_params, **mask_kwargs)
    ellip_list = getattr(data, 'get_{0}_ellipse_list'.format(kind))()
    for e in ellip_list:
        ax.add_artist(e)
    return ax


def plot_original(data, grid, sub_grid_ratio=[0.95, 0.05]):
    gs_inner = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=sub_grid_ratio, subplot_spec=grid, wspace=0.01)
    ax = make_ax(gs_inner[0], data.wcs, color_grid='C7')
    ax.add_patch(data.get_hexagon(correct_hex=True, edgecolor='C4'))
    ax.imshow(data.image)
    return ax


def plot_spectra(data, grid):
    ax = plt.subplot(grid)
    mean_spiral, mean_bar, mean_center = get_stacked_spectra(data)
    if mean_center is not None:
        ax.plot(mean_center.wavelength, mean_center.flux, color='C2', label='center', alpha=0.8)
    if mean_bar is not None:
        ax.plot(mean_bar.wavelength, mean_bar.flux, color='C1', label='bar', alpha=0.8)
    if mean_spiral is not None:
        ax.plot(mean_spiral.wavelength, mean_spiral.flux, color='C0', label='spiral', alpha=0.8)
    ax.set_xlabel(r'Wavelength $[\rm\AA]$')
    ax.set_ylabel(r'Flux $[\rm 10^{-17}\,erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$')
    ax.set_ylim([0, None])
    ax.set_title('Stacked spectra')
    ax.legend(loc=1)
    return ax


def gz3d_plot(filename, fdx, output_folder, output_file, spaxel_grid=False):
    if (os.path.isfile('{0}/no_mask/{1}.png'.format(output_folder, output_file))) or (os.path.isfile('{0}/mask/{1}.png'.format(output_folder, output_file))):
        # don't re-make plots if they already exist
        return 'already_exists'
    data = gz3d_fits(filename)
    fig_width = 16.4
    fig_height = 13.5
    sub_grid_ratio = [0.95, 0.05]
    gs = gridspec.GridSpec(3, 3)
    gs.update(left=0.05, right=0.94, bottom=0.05, top=0.94, wspace=0.5, hspace=0.3)
    fig = plt.figure(fdx, figsize=(fig_width, fig_height))
    # plot original image
    ax = plot_original(data, gs[0, 0], sub_grid_ratio=sub_grid_ratio)
    ax.set_title(' \n {0}'.format(data.metadata['MANGAID'][0]))
    # plot center image
    ax = plot_ellipse(data, 'center', gs[0, 1], sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
    ax.set_title('Center(s) \n {0} Classifications'.format(data.num_center_star_classifications))
    # plot star image
    ax = plot_ellipse(data, 'star', gs[0, 2], sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
    ax.set_title('Star(s) \n {0} Classifications'.format(data.num_center_star_classifications))
    # plot spiral image
    ax = plot_mask(data, 'spiral', gs[1, 0], sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
    ax.set_title('Spiral arms \n {0} Classifications'.format(data.num_spiral_classifications))
    # plot bar image
    ax = plot_mask(data, 'bar', gs[1, 1], sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
    ax.set_title('Bar \n {0} Classifications'.format(data.num_bar_classifications))
    # save figure and close data
    subfolder = 'no_mask'
    if (data.bar_mask.sum() > 0) or (data.spiral_mask.sum() > 0):
        ax = plot_spectra(data, gs[2, :])
        subfolder = 'mask'
    output_file = '{0}/{1}/{2}.png'.format(output_folder, subfolder, output_file)
    fig.savefig(output_file)
    data.close()
    plt.close(fig)
    return 'done'


def make_plots(fits_location, output_folder):
    file_list = os.listdir(fits_location)
    widgets = ['Plot: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
    pbar = pb.ProgressBar(widgets=widgets, maxval=len(file_list))
    pbar.start()
    for fdx, filename in enumerate(file_list):
        if '.fits.gz' in filename:
            output_file = filename.split('.')[0]
            status = gz3d_plot(os.path.join(fits_location, filename), fdx, output_folder, output_file)
        pbar.update(fdx + 1)
    pbar.finish()


if __name__ == '__main__':
    import warnings
    warnings.filterwarnings('ignore')
    make_plots('/Volumes/Work/GZ3D/MPL5_fits', '/Volumes/Work/GZ3D/MPL5_plots_v3')

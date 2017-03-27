from astropy.io import fits
from astropy.wcs import WCS
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as patches
import mpl_style
import os
import progressbar as pb
import scipy.linalg as sl

plt.style.use(mpl_style.style1)
LUT = {7: 3, 19: 5, 37: 7, 61: 9, 91: 11, 127: 13}


def get_hexagon(ra, dec, ifu_size, wcs):
    diameter = 2.0 / 0.099
    c = wcs.wcs_world2pix(np.array([[ra, dec]]), 1)
    r = LUT[ifu_size] * diameter * np.sqrt(3) / 4
    return patches.RegularPolygon(c[0], 6, r, fill=False, orientation=np.deg2rad(30), edgecolor='magenta', linewidth=0.8)


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


def cov_to_ellipse(cov, pos, nstd=1, **kwargs):
    eigvec, eigval, V = sl.svd(cov, full_matrices=False)
    # the angle the first eigenvector makes with the x-axis
    theta = np.degrees(np.arctan2(eigvec[1, 0], eigvec[0, 0]))
    # full width and height of ellipse, not radius
    # the eigenvalues are the variance along the eigenvectors
    width, height = 2 * nstd * np.sqrt(eigval)
    return patches.Ellipse(xy=pos, width=width, height=height, angle=theta, **kwargs)


def plot_mask(hdu, wcs, grid, title, hexparams, cmap=plt.cm.Purples, sub_grid_ratio=[1, 0.2]):
    gs_inner = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=sub_grid_ratio, subplot_spec=grid)
    ax = plt.subplot(gs_inner[0], projection=wcs)
    ax.set_title(title)
    set_up_axes(ax, color_grid='#1f77b4', color_tick='black')
    ax.add_patch(get_hexagon(*hexparams))
    vmax = hdu.data.max()
    if vmax != 0:
        im = ax.imshow(hdu.data, cmap=cmap, vmin=0, vmax=vmax)
        ax_bar = plt.subplot(gs_inner[1])
        cb = plt.colorbar(im, cax=ax_bar)
        cb.set_label('Count')
    else:
        im = ax.imshow(hdu.data, cmap=cmap)
    return ax


def plot_ellipse(hdu, *mask_params, **mask_kwargs):
    ax = plot_mask(*mask_params, **mask_kwargs)
    for idx in range(len(hdu.data)):
        pos = np.array([hdu.data['x'][idx], hdu.data['y'][idx]])
        cov = np.array([[hdu.data['var_x'][idx], hdu.data['var_x_y'][idx]], [hdu.data['var_x_y'][idx], hdu.data['var_y'][idx]]])
        ellip = cov_to_ellipse(cov, pos, nstd=2, edgecolor='k', facecolor='none', lw=1)
        ax.add_artist(ellip)


def plot_original(hdu, wcs, grid, title, sub_grid_ratio=[1, 0.2]):
    gs_inner = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=sub_grid_ratio, subplot_spec=grid)
    ax = plt.subplot(gs_inner[0], projection=wcs)
    ax.set_title(title)
    set_up_axes(ax, color_grid='#1f77b4')
    ax.imshow(hdu.data)


def gz3d_plot(filename, fdx, output_file):
    hdulist = fits.open(filename)
    fig_width = 16.4
    fig_height = 9
    sub_grid_ratio = [0.95, 0.05]
    gs = gridspec.GridSpec(2, 3)
    gs.update(left=0.05, right=0.94, bottom=0.05, top=0.94, wspace=0.5, hspace=0.3)
    fig = plt.figure(fdx, figsize=(fig_width, fig_height))
    wcs = WCS(hdulist[1].header)
    hexparams = (hdulist[5].data['ra'][0], hdulist[5].data['dec'][0], int(hdulist[5].data['IFUDESIGNSIZE'][0]), wcs)
    # plot original image
    title = ' \n {0}'.format(hdulist[5].data['MANGAID'][0])
    plot_original(hdulist[0], wcs, gs[0, 0], title, sub_grid_ratio=sub_grid_ratio)
    # plot center image
    title = 'Center(s) \n {0} Classifications'.format(len(hdulist[8].data))
    plot_ellipse(hdulist[6], hdulist[1], wcs, gs[0, 1], title, hexparams, sub_grid_ratio=sub_grid_ratio)
    # plot star image
    title = 'Star(s) \n {0} Classifications'.format(len(hdulist[8].data))
    plot_ellipse(hdulist[7], hdulist[2], wcs, gs[0, 2], title, hexparams, sub_grid_ratio=sub_grid_ratio)
    # plot spiral image
    title = 'Spiral arms \n {0} Classifications'.format(len(hdulist[9].data))
    plot_mask(hdulist[3], wcs, gs[1, 0], title, hexparams, sub_grid_ratio=sub_grid_ratio)
    # plot bar image
    title = 'Bar \n {0} Classifications'.format(len(hdulist[10].data))
    plot_mask(hdulist[4], wcs, gs[1, 1], title, hexparams, sub_grid_ratio=sub_grid_ratio)
    fig.savefig('{0}.png'.format(output_file))
    plt.close(fig)


def make_plots(fits_location, output_folder):
    file_list = os.listdir(fits_location)
    widgets = ['Plot: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
    pbar = pb.ProgressBar(widgets=widgets, maxval=len(file_list))
    pbar.start()
    for fdx, filename in enumerate(file_list):
        if '.fits.gz' in filename:
            output_file = '{0}/{1}'.format(output_folder, filename.split('.')[0])
            gz3d_plot(os.path.join(fits_location, filename), fdx, output_file)
        pbar.update(fdx + 1)
    pbar.finish()


if __name__ == '__main__':
    make_plots('/Volumes/Work/GZ3D/MPL5_fits', '/Volumes/Work/GZ3D/MPL5_plots')

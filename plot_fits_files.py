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
from scipy.interpolate import RectBivariateSpline

plt.style.use(mpl_style.style1)
LUT = {7: 3, 19: 5, 37: 7, 61: 9, 91: 11, 127: 13}
spaxel_grid = {7: 23, 19: 33, 37: 43, 61: 53, 91: 63, 127: 75}


def get_hexagon(ra, dec, ifu_size, wcs, correct_hex=False, edgecolor='magenta'):
    # the spacing should be ~0.5 arcsec not 0, and it should not be rotated by np.sqrt(3) / 2
    if correct_hex:
        # each hex has a total diameter of 2.5 arcsec on the sky (only 2 of it is a fiber)
        diameter = 2.5 / 0.099
        # the radius for mpl is from the center to each vertex, not center to side
        r = LUT[ifu_size] * diameter / 2
    else:
        # this was me being wrong about all the hexagon params
        # these hexagons are about 0.7 times too small (2 * np.sqrt(3) / 5 to be exact)
        diameter = 2.0 / 0.099
        r = LUT[ifu_size] * diameter * np.sqrt(3) / 4
    c = wcs.wcs_world2pix(np.array([[ra, dec]]), 1)
    return patches.RegularPolygon(c[0], 6, r, fill=False, orientation=np.deg2rad(30), edgecolor=edgecolor, linewidth=0.8)


def get_spaxel_grid(ra, dec, ifu_size, wcs, ax):
    # grid size in arcsec
    one_grid = 0.5 / 0.099
    c = wcs.wcs_world2pix(np.array([[ra, dec]]), 1)
    grid = np.arange(spaxel_grid[ifu_size] + 1) * one_grid
    grid_x = grid - np.median(grid) + c[0, 0]
    grid_y = grid - np.median(grid) + c[0, 1]
    v_line_x = np.array(zip(grid_x, grid_x)).T
    v_line_y = np.array([(grid_y[0], grid_y[-1])]).T
    h_line_x = np.array([(grid_x[0], grid_x[-1])]).T
    h_line_y = np.array(zip(grid_y, grid_y)).T
    ax.plot(v_line_x, v_line_y, color='C7', lw=0.5, alpha=0.3)
    ax.plot(h_line_x, h_line_y, color='C7', lw=0.5, alpha=0.3)


def spaxel_mask(mask, ra, dec, ifu_size, wcs):
    # assues a 0.5 arcsec grid centered on the ifu's ra and dec
    # use a Bivariate spline approximation to resample maks to the spaxel grid
    xx = np.arange(mask.shape[1])
    yy = np.arange(mask.shape[0])
    s = RectBivariateSpline(xx, yy, mask)
    one_grid = 0.5 / 0.099
    c = wcs.wcs_world2pix(np.array([[ra, dec]]), 1)
    grid = np.arange(spaxel_grid[ifu_size]) * one_grid
    grid_x = grid - np.median(grid) + c[0, 0]
    grid_y = grid - np.median(grid) + c[0, 1]
    # flip the output mask so the origin is the lower left of the image
    s_mask = np.flipud(s(grid_x, grid_y))
    # zero out small values
    s_mask[s_mask < 0.5] = 0
    return s_mask


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


def plot_mask(hdu, wcs, grid, title, hexparams, cmap=plt.cm.Purples, sub_grid_ratio=[1, 0.2], spaxel_grid=False):
    gs_inner = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=sub_grid_ratio, subplot_spec=grid)
    ax = plt.subplot(gs_inner[0], projection=wcs)
    ax.set_title(title)
    set_up_axes(ax, color_grid='C0', color_tick='black')
    ax.add_patch(get_hexagon(*hexparams))
    ax.add_patch(get_hexagon(*hexparams, correct_hex=True, edgecolor='C4'))
    if spaxel_grid:
        get_spaxel_grid(hexparams[0], hexparams[1], hexparams[2], hexparams[3], ax)
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


def plot_original(hdu, wcs, grid, title, hexparams, sub_grid_ratio=[1, 0.2]):
    gs_inner = gridspec.GridSpecFromSubplotSpec(1, 2, width_ratios=sub_grid_ratio, subplot_spec=grid)
    ax = plt.subplot(gs_inner[0], projection=wcs)
    ax.set_title(title)
    set_up_axes(ax, color_grid='C0')
    ax.add_patch(get_hexagon(*hexparams, correct_hex=True, edgecolor='C4'))
    ax.imshow(hdu.data)


def gz3d_plot(filename, fdx, output_file, spaxel_grid=False):
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
    plot_original(hdulist[0], wcs, gs[0, 0], title, hexparams, sub_grid_ratio=sub_grid_ratio)
    # plot center image
    title = 'Center(s) \n {0} Classifications'.format(len(hdulist[8].data))
    plot_ellipse(hdulist[6], hdulist[1], wcs, gs[0, 1], title, hexparams, sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
    # plot star image
    title = 'Star(s) \n {0} Classifications'.format(len(hdulist[8].data))
    plot_ellipse(hdulist[7], hdulist[2], wcs, gs[0, 2], title, hexparams, sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
    # plot spiral image
    title = 'Spiral arms \n {0} Classifications'.format(len(hdulist[9].data))
    plot_mask(hdulist[3], wcs, gs[1, 0], title, hexparams, sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
    # plot bar image
    title = 'Bar \n {0} Classifications'.format(len(hdulist[10].data))
    plot_mask(hdulist[4], wcs, gs[1, 1], title, hexparams, sub_grid_ratio=sub_grid_ratio, spaxel_grid=spaxel_grid)
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

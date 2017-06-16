import matplotlib as mpl
import numpy as np
import astropy.wcs as wcs


def plot_alpha_scatter(x, y, mask, color, ax, snr=None, sf_mask=None, value=True, **kwargs):
    idx = mask > 0
    if value:
        idx = idx & (y.value > 0)
    if (value) and (snr is not None):
        idx = idx & (y.snr > snr)
    if sf_mask is not None:
        idx = idx & sf_mask
    c = mpl.colors.to_rgb(color)
    c_a = np.array([c + (i, ) for i in mask[idx] / 15])
    c_a[c_a > 1] = 1
    if value:
        return ax.scatter(x[idx], y.value[idx], c=c_a, edgecolor=c_a, **kwargs)
    else:
        return ax.scatter(x[idx], y[idx], c=c_a, edgecolor=c_a, **kwargs)


def plot_by_r(gz3d, ax, key='specindex_dn4000', ylabel=r'$D_{n}4000$', snr=3, sf_only=False, s=15):
    title = 'S/N > {0}'.format(snr)
    r = gz3d.maps['spx_ellcoo_elliptical_radius'].value
    r_50 = gz3d.maps.nsa['elpetro_th50_r']
    line = gz3d.maps[key]
    sf_mask = None
    if sf_only:
        sf_mask = gz3d.sf_mask
        title += ', star forming only'
    # plot scatter points
    plot_alpha_scatter(r / r_50, line, gz3d.spiral_mask_spaxel, 'C0', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(r / r_50, line, gz3d.bar_mask_spaxel, 'C1', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(r / r_50, line, gz3d.star_mask_spaxel, 'C3', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(r / r_50, line, gz3d.center_mask_spaxel, 'C2', ax, s=s, snr=snr, sf_mask=sf_mask)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r'R / R$_{50}$')


def plot_by_theta(gz3d, ax, key='specindex_dn4000', ylabel=r'$D_{n}4000$', snr=3, sf_only=False, s=15):
    title = 'S/N > {0}'.format(snr)
    theta = gz3d.maps['spx_ellcoo_elliptical_azimuth'].value
    line = gz3d.maps[key]
    sf_mask = None
    if sf_only:
        sf_mask = gz3d.sf_mask
        title += ', star forming only'
    # plot scatter points
    plot_alpha_scatter(theta, line, gz3d.spiral_mask_spaxel, 'C0', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(theta, line, gz3d.bar_mask_spaxel, 'C1', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(theta, line, gz3d.star_mask_spaxel, 'C3', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(theta, line, gz3d.center_mask_spaxel, 'C2', ax, s=s, snr=snr, sf_mask=sf_mask)
    ax.set_xticks([0, 90, 180, 270, 360])
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r'$\theta$')


def zero_theta_line(gz3d):
    phi = gz3d.maps.nsa['elpetro_phi']
    ra = gz3d.cube.ra
    dec = gz3d.cube.dec
    map_wcs = wcs.WCS(gz3d.maps['spx_ellcoo_elliptical_azimuth'].header, naxis=2)
    # get the center of the image
    cx, cy = map_wcs.wcs_world2pix(ra, dec, 0)
    # get the max radius
    r = np.sqrt(cx**2 + cy**2)
    # get the end of the line
    x = r * np.sin(np.deg2rad(-phi)) + cx
    y = r * np.cos(np.deg2rad(-phi)) + cy
    # world coords
    ra, dec = map_wcs.wcs_pix2world([cx, x], [cy, y], 0)
    # image coords
    return gz3d.wcs.wcs_world2pix(ra, dec, 0)


if __name__ == '__main__':
    from gz3d_fits import gz3d_fits
    import matplotlib.pyplot as plt
    file_name = '/Volumes/Work/GZ3D/MPL5_fits/1-167242_127_5679242.fits.gz'
    gz3d = gz3d_fits(file_name)
    gz3d.get_cube(maps=True)
    gz3d.make_all_spaxel_masks()
    gz3d.get_bpt()
    plt.figure(1)
    ax1 = plt.gca()
    plot_by_r(gz3d, ax1, key='specindex_dn4000', ylabel=r'$D_{n}4000$')
    plt.figure(2)
    ax2 = plt.gca()
    plot_by_r(gz3d, ax2, key='emline_sew_ha_6564', ylabel=r'EW(H$\alpha$)', sf_only=True)
    plt.figure(3)
    ax1 = plt.gca()
    plot_by_theta(gz3d, ax1, key='specindex_dn4000', ylabel=r'$D_{n}4000$')
    plt.figure(4)
    ax2 = plt.gca()
    plot_by_theta(gz3d, ax2, key='emline_sew_ha_6564', ylabel=r'EW(H$\alpha$)', sf_only=True)
    plt.show()

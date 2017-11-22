import matplotlib as mpl
import numpy as np
import astropy.wcs as wcs
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt


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
    r_plot = r / r_50
    line = gz3d.maps[key]
    # other spaxel masks
    odx = (gz3d.other_mask_spaxel > 0) & (line.value > 0)
    if snr is not None:
        odx = odx & (line.snr > snr)
    sf_mask = None
    if sf_only:
        sf_mask = gz3d.sf_mask
        title += ', star forming only'
        odx = odx & sf_mask
    # plot scatter points
    ax.scatter(r_plot[odx], line.value[odx], c='#c5c5c5', edgecolor='#c5c5c5', s=s)
    plot_alpha_scatter(r_plot, line, gz3d.spiral_mask_spaxel, 'C0', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(r_plot, line, gz3d.bar_mask_spaxel, 'C1', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(r_plot, line, gz3d.star_mask_spaxel, 'C3', ax, s=s, snr=snr, sf_mask=sf_mask)
    plot_alpha_scatter(r_plot, line, gz3d.center_mask_spaxel, 'C2', ax, s=s, snr=snr, sf_mask=sf_mask)
    ax.set_title(title)
    ax.set_ylabel(ylabel)
    ax.set_xlabel(r'R / R$_{50}$')


def plot_by_theta(gz3d, ax, key='specindex_dn4000', ylabel=r'$D_{n}4000$', snr=3, sf_only=False, s=15):
    title = 'S/N > {0}'.format(snr)
    theta = gz3d.maps['spx_ellcoo_elliptical_azimuth'].value
    line = gz3d.maps[key]
    # other spaxel masks
    odx = (gz3d.other_mask_spaxel > 0) & (line.value > 0)
    if snr is not None:
        odx = odx & (line.snr > snr)
    sf_mask = None
    if sf_only:
        sf_mask = gz3d.sf_mask
        title += ', star forming only'
        odx = odx & sf_mask
    # plot scatter points
    ax.scatter(theta[odx], line.value[odx], c='#c5c5c5', edgecolor='#c5c5c5', s=s)
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
    map_wcs = wcs.WCS(gz3d.maps['spx_ellcoo_elliptical_azimuth'].header, naxis=2)
    # get the center of the image
    # cx, cy = map_wcs.wcs_world2pix(ra, dec, 0)
    r_map = gz3d.maps['spx_ellcoo_elliptical_radius'].value
    cy, cx = map(np.mean, np.nonzero(np.isclose(r_map, 0, atol=1)))
    # get the max radius
    r = np.sqrt(cx**2 + cy**2)
    # get the end of the line
    x = r * np.sin(np.deg2rad(-phi)) + cx
    y = r * np.cos(np.deg2rad(-phi)) + cy
    # world coords
    ra_line, dec_line = map_wcs.wcs_pix2world([cx, x], [cy, y], 0)
    # image coords
    return gz3d.wcs.wcs_world2pix(ra_line, dec_line, 0)


def r_theta_plot(gz3d, fdx=0, **kwargs):
    fig_width = 12
    fig_height = 4.5
    gs = gridspec.GridSpec(1, 2)
    gs.update(left=0.05, right=0.94, bottom=0.05, top=0.94, wspace=0.4, hspace=0.3)
    fig = plt.figure(fdx, figsize=(fig_width, fig_height))
    ax_00 = plt.subplot(gs[0, 0])
    plot_by_r(gz3d, ax_00, **kwargs)
    ax_01 = plt.subplot(gs[0, 1])
    plot_by_theta(gz3d, ax_01, **kwargs)
    return fig


if __name__ == '__main__':
    from gz3d_fits import gz3d_fits
    import mpl_style
    import os
    plt.style.use(mpl_style.style1)
    '''
    file_name = '/Volumes/Work/GZ3D/MPL5_fits/1-167242_127_5679242.fits.gz'
    gz3d = gz3d_fits(file_name)
    gz3d.get_cube(maps=True)
    gz3d.make_all_spaxel_masks()
    gz3d.get_bpt()
    plt.figure(1)
    ax1 = plt.gca()
    plot_by_r(gz3d, ax1, key='specindex_dn4000', ylabel=r'$D_{n}4000$', s=8)
    plt.figure(2)
    ax2 = plt.gca()
    plot_by_r(gz3d, ax2, key='emline_sew_ha_6564', ylabel=r'EW(H$\alpha$)', sf_only=True, s=8)
    plt.figure(3)
    ax1 = plt.gca()
    plot_by_theta(gz3d, ax1, key='specindex_dn4000', ylabel=r'$D_{n}4000$', s=8)
    plt.figure(4)
    ax2 = plt.gca()
    plot_by_theta(gz3d, ax2, key='emline_sew_ha_6564', ylabel=r'EW(H$\alpha$)', sf_only=True, s=8)
    plt.show()
    '''
    filepath = '/Volumes/Work/GZ3D/MPL5_fits'
    output_folder_dn = '/Users/coleman/Desktop/plots_for_talk/Dn_4000'
    output_folder_ha = '/Users/coleman/Desktop/plots_for_talk/H_alpha'
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
    fdx = 0
    for mid in id_list:
        output_name_dn = os.path.join(output_folder_dn, mid) + '_dn.png'
        output_name_ha = os.path.join(output_folder_ha, mid) + '_ha.png'
        input_name = os.path.join(filepath, mid) + '.fits.gz'
        gz3d = gz3d_fits(input_name)
        gz3d.get_bpt()
        gz3d.make_all_spaxel_masks()
        fig_dn = r_theta_plot(gz3d, fdx=fdx, key='specindex_dn4000', ylabel=r'$D_{n}4000$', s=8)
        fig_dn.savefig(output_name_dn)
        plt.close(fig_dn)
        fig_ha = r_theta_plot(gz3d, fdx=fdx, key='emline_sew_ha_6564', ylabel=r'EW(H$\alpha$)', sf_only=True, s=8)
        fig_ha.savefig(output_name_ha)
        plt.close(fig_ha)

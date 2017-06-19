from gz3d_fits import gz3d_fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import mpl_style
import plot_fits_files
from alpha_overlay_maps import plot_masks
from plot_fits_with_bpt import plot_bpt, plot_bpt_boundry, plot_spectra, make_spectra_ax
from plot_by_r import plot_by_r, plot_by_theta

plt.style.use(mpl_style.style1)


def diagnostic_plots(filename, output_name=None, fdx=1, oi_sf=False):
    gz3d = gz3d_fits(filename)
    gz3d.get_bpt(oi_sf=oi_sf)
    gz3d.make_all_spaxel_masks()
    fig_width = 21.9
    fig_height = 13.5
    gs = gridspec.GridSpec(3, 4)
    gs.update(left=0.05, right=0.94, bottom=0.05, top=0.94, wspace=0.4, hspace=0.3)
    fig = plt.figure(fdx, figsize=(fig_width, fig_height))
    # plot original image
    ax_00 = plot_fits_files.plot_original(gz3d, gs[0, 0], sub_grid_ratio=[0.9, 0.1])
    ax_00.set_title('{0}'.format(gz3d.metadata['MANGAID'][0]))
    # plot masks
    ax_10 = plot_masks(gz3d, gs[1, 0], sub_grid_ratio=[0.9, 0.1])
    # plot dn4000
    ax_01 = plt.subplot(gs[0, 1])
    plot_by_r(gz3d, ax_01, key='specindex_dn4000', ylabel=r'$D_{n}4000$', s=8)
    ax_11 = plt.subplot(gs[1, 1])
    plot_by_theta(gz3d, ax_11, key='specindex_dn4000', ylabel=r'$D_{n}4000$', s=8)
    # plot EW(Ha)
    ax_01 = plt.subplot(gs[0, 2])
    plot_by_r(gz3d, ax_01, key='emline_sew_ha_6564', ylabel=r'EW(H$\alpha$)', sf_only=True, s=8)
    ax_11 = plt.subplot(gs[1, 2])
    plot_by_theta(gz3d, ax_11, key='emline_sew_ha_6564', ylabel=r'EW(H$\alpha$)', sf_only=True, s=8)
    # plot bpt
    axs = [plt.subplot(gs[0, 3]), plt.subplot(gs[1, 3]), plt.subplot(gs[2, 3])]
    plot_bpt(axs, gz3d, s=3)
    plot_bpt_boundry(axs)
    # plot spectra
    gz3d.get_mean_spectra(inv=True)
    ax_20 = plt.subplot(gs[2, :-1])
    plot_spectra(ax_20, gz3d, 'center', label=r'$\Delta$Center', color='C2', alpha=0.8)
    plot_spectra(ax_20, gz3d, 'bar', label=r'$\Delta$Bar', color='C1', alpha=0.8)
    plot_spectra(ax_20, gz3d, 'spiral', label=r'$\Delta$Spiral', color='C0', alpha=0.8)
    make_spectra_ax(ax_20)
    if output_name is not None:
        fig.savefig('{0}.png'.format(output_name))
    gz3d.close()
    plt.close(fig)


if __name__ == '__main__':
    file_name = '/Volumes/Work/GZ3D/MPL5_fits/1-167242_127_5679242.fits.gz'
    diagnostic_plots(file_name)
    plt.show()

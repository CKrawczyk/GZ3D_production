from gz3d_fits import gz3d_fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import to_rgb
import marvin.utils.dap.bpt as bpt
import mpl_style
import plot_fits_files
import numpy as np
from plot_by_r import plot_alpha_scatter

plt.style.use(mpl_style.style1)


def plot_bpt_boundry(axs):
    # Plots the classification boundary lines
    xx_sf_nii = np.linspace(-1.281, 0.045, int(1e4))
    xx_sf_sii = np.linspace(-2, 0.315, int(1e4))
    xx_sf_oi = np.linspace(-2.5, -0.7, int(1e4))
    xx_comp_nii = np.linspace(-2, 0.4, int(1e4))
    xx_agn_sii = np.array([-0.308, 1.0])
    xx_agn_oi = np.array([-1.12, 0.5])
    # xlims = [ax.get_xlim() for ax in axs]
    # ylims = [ax.get_ylim() for ax in axs]
    xlims = [(-2, 0.5), (-1.5, 0.5), (-2.5, 0.0)]
    ylims = [(-1.5, 1.3), (-1.5, 1.3), (-1.5, 1.3)]
    xlabels = [r'log([NII]/H$\alpha$)', r'log([SII]/H$\alpha$)', r'log([OI]/H$\alpha$)']
    axs[0].plot(xx_sf_nii, bpt.kewley_sf_nii(xx_sf_nii), 'k--', zorder=90)
    axs[0].plot(xx_comp_nii, bpt.kewley_comp_nii(xx_comp_nii), 'k-', zorder=90)
    axs[1].plot(xx_sf_sii, bpt.kewley_sf_sii(xx_sf_sii), 'k-', zorder=90)
    axs[1].plot(xx_agn_sii, bpt.kewley_agn_sii(xx_agn_sii), 'k-', zorder=80)
    axs[2].plot(xx_sf_oi, bpt.kewley_sf_oi(xx_sf_oi), 'k-', zorder=90)
    axs[2].plot(xx_agn_oi, bpt.kewley_agn_oi(xx_agn_oi), 'k-', zorder=80)
    for ax, xlim, ylim, xlabel in zip(axs, xlims, ylims, xlabels):
        ax.set_ylabel(r'log([OIII]/H$\beta$)')
        ax.set_xlabel(xlabel)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
    axs[0].legend(loc=3)


def plot_bpt_alpha(gz3d, ax, bpt_kind, **kwargs):
    y = gz3d.log_oiii_hb
    x = getattr(gz3d, bpt_kind)
    mdx = ~(y.mask | x.mask)
    s = kwargs.pop('s', 8)
    odx = mdx & (gz3d.other_mask_spaxel > 0)
    ax.scatter(x[odx], y[odx], c='#c5c5c5', edgecolor='#c5c5c5', s=s, label='Other', **kwargs)
    plot_alpha_scatter(x, y, gz3d.spiral_mask_spaxel, 'C0', ax, s=s, sf_mask=mdx, snr=None, value=False, **kwargs)
    plot_alpha_scatter(x, y, gz3d.bar_mask_spaxel, 'C1', ax, s=s, sf_mask=mdx, snr=None, value=False, **kwargs)
    plot_alpha_scatter(x, y, gz3d.star_mask_spaxel, 'C3', ax, s=s, sf_mask=mdx, snr=None, value=False, **kwargs)
    plot_alpha_scatter(x, y, gz3d.center_mask_spaxel, 'C2', ax, s=s, sf_mask=mdx, snr=None, value=False, **kwargs)


def plot_bpt(axs, gz3d, **kwargs):
    for idx, bpt_name in enumerate(['log_nii_ha', 'log_sii_ha', 'log_oi_ha']):
        plot_bpt_alpha(gz3d, axs[idx], bpt_name, **kwargs)


def screen_color(c1, c2):
    c1 = np.array(c1)
    c2 = np.array(c2)
    return 1 - ((1 - c1) * (1 - c2))


def plot_bpt_dis(axs, gz3d, mask_name, color='C0', **kwargs):
    if gz3d.cube != 'no_data':
        for idx, bpt_name in enumerate(['log_nii_ha', 'log_sii_ha', 'log_oi_ha']):
            x, y, dis = gz3d.bpt_in_mask(mask_name, bpt_name)
            if len(dis) > 0:
                rgb = to_rgb(color)
                colors = np.array([screen_color(rgb, d) for d in dis])
                axs[idx].scatter(x, y, color=colors, edgecolor=colors, **kwargs)


def plot_spectra(ax, gz3d, spectrum_name, **kwargs):
    spectrum = getattr(gz3d, 'mean_{0}'.format(spectrum_name))
    spectrum_inv = getattr(gz3d, 'mean_not_{0}'.format(spectrum_name))
    if spectrum is not None:
        spectrum_dif = spectrum - spectrum_inv
        ax.plot(spectrum_dif.wavelength, spectrum_dif.flux, **kwargs)


def make_spectra_ax(ax):
    ax.set_xlabel(r'Wavelength $[\rm\AA]$')
    ax.set_ylabel(r'$\Delta$Flux $[\rm 10^{-17}\,erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$')
    # ax.set_ylim([0, None])
    ax.set_title('Difference spectra')
    ax.legend(loc=1)


if __name__ == '__main__':
    from astropy.table import Table
    import os
    import progressbar as pb
    filepath = '/Volumes/Work/GZ3D/MPL5_fits'
    filenames = Table.read('bar_or_spiral.txt', format='ascii.csv')
    output_folder = '/Volumes/Work/GZ3D/MPL5_plots_bpt'
    widgets = ['Plot: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
    pbar = pb.ProgressBar(widgets=widgets, maxval=len(filenames))
    pbar.start()
    for fdx, filename in enumerate(filenames['filename']):
        output_name = filename.split('.')[0]
        # skip existing files
        if os.path.isfile('{0}/{1}.png'.format(output_folder, output_name)):
            pbar.update(fdx + 1)
            continue
        gz3d = gz3d_fits('{0}/{1}'.format(filepath, filename))
        gz3d.get_bpt()
        gz3d.get_mean_spectra(inv=True)
        fig_width = 21.9
        fig_height = 13.5
        gs = gridspec.GridSpec(3, 4)
        gs.update(left=0.05, right=0.94, bottom=0.05, top=0.94, wspace=0.5, hspace=0.3)
        fig = plt.figure(fdx, figsize=(fig_width, fig_height))
        # plot original image
        ax_00 = plot_fits_files.plot_original(gz3d, gs[0, 0])
        ax_00.set_title('{0}'.format(gz3d.metadata['MANGAID'][0]))
        # plot center image
        ax_01 = plot_fits_files.plot_ellipse(gz3d, 'center', gs[0, 1])
        ax_01.set_title('Center(s) \n {0} Classifications'.format(gz3d.num_center_star_classifications))
        # plot star image
        ax_02 = plot_fits_files.plot_ellipse(gz3d, 'star', gs[0, 2])
        ax_02.set_title('Star(s) \n {0} Classifications'.format(gz3d.num_center_star_classifications))
        # plot spiral image
        ax_10 = plot_fits_files.plot_mask(gz3d, 'spiral', gs[1, 0])
        ax_10.set_title('Spiral arms \n {0} Classifications'.format(gz3d.num_spiral_classifications))
        # plot bar image
        ax_11 = plot_fits_files.plot_mask(gz3d, 'bar', gs[1, 1])
        ax_11.set_title('Bar \n {0} Classifications'.format(gz3d.num_bar_classifications))
        # plot spectra
        ax_20 = plt.subplot(gs[2, :-1])
        plot_spectra(ax_20, gz3d, 'center', label=r'$\Delta$Center', color='C2', alpha=0.8)
        plot_spectra(ax_20, gz3d, 'bar', label=r'$\Delta$Bar', color='C1', alpha=0.8)
        plot_spectra(ax_20, gz3d, 'spiral', label=r'$\Delta$Spiral', color='C0', alpha=0.8)
        make_spectra_ax(ax_20)
        # plot bpt
        axs = [plt.subplot(gs[0, 3]), plt.subplot(gs[1, 3]), plt.subplot(gs[2, 3])]
        plot_bpt(axs, gz3d, 'other_mask_spaxel', color='C7', label='Other', s=3)
        plot_bpt(axs, gz3d, 'spiral_mask_spaxel', color='C0', label='Spiral', s=3)
        plot_bpt(axs, gz3d, 'bar_mask_spaxel', color='C1', label='Bar', s=3)
        plot_bpt(axs, gz3d, 'center_mask_spaxel', color='C2', label='Center', s=3)
        plot_bpt_boundry(axs)
        # save plot
        fig.savefig('{0}/{1}.png'.format(output_folder, output_name))
        # close files and figure
        gz3d.close()
        plt.close(fig)
        pbar.update(fdx + 1)
    pbar.finish()

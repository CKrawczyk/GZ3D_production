from marvin import config
from marvin_subclass import CubeFast
from gz3d_fits import gz3d_fits
import matplotlib.pyplot as plt
import numpy as np
import sys
import traceback


class Suppressor(object):
    # A class to mute the output of a function
    def __enter__(self):
        self.stdout = sys.stdout
        sys.stdout = self

    def __exit__(self, type, value, traceback):
        sys.stdout = self.stdout

    def write(self, x):
        pass


# set marvin to local mode
# downloading the cubes gives faster access to the spectra
config.mode = 'local'
config.download = True


def mean_spectrum(spectra, weights=None):
    if weights is None:
        weights = np.ones(len(spectra))
    inside_ifu = np.array([s.inside_ifu() for s in spectra])
    stack = (spectra[inside_ifu] * weights[inside_ifu]).sum()
    return stack / weights[inside_ifu].sum()


def stack_spectra(cube, mask):
    yy, xx = np.where(mask > 0)
    if len(yy) > 0:
        mask_value = mask[yy, xx]
        spaxels = cube[yy, xx]
        spectra = np.array([s.spectrum for s in spaxels])
        if len(spectra) == 1:
            return spectra[0]
        else:
            return mean_spectrum(spectra, weights=mask_value)
    else:
        return None


def inv_mask(mask):
    return mask.max() - mask


def get_stacked_spectra(gz3d, cube=None, return_inv=False):
    manga_id = gz3d.metadata['MANGAID'][0]
    output = []
    try:
        if cube is not None:
            c = cube
        else:
            with Suppressor():
                # mute marvin so there are no `SDSS_ACCESS` prints to stdout
                c = CubeFast(mangaid=manga_id)
        gz3d.make_all_spaxel_masks(grid_size=c.data['FLUX'].data.shape[1:])
        mean_bar = stack_spectra(c, gz3d.bar_mask_spaxel)
        mean_spiral = stack_spectra(c, gz3d.spiral_mask_spaxel)
        mean_center = stack_spectra(c, gz3d.center_mask_spaxel)
        output += [mean_spiral, mean_bar, mean_center]
        if return_inv:
            mean_not_bar = stack_spectra(c, inv_mask(gz3d.bar_mask_spaxel))
            mean_not_spiral = stack_spectra(c, inv_mask(gz3d.spiral_mask_spaxel))
            mean_not_center = stack_spectra(c, inv_mask(gz3d.center_mask_spaxel))
            output += [mean_not_spiral, mean_not_bar, mean_not_center]
        if cube is None:
            c.data.close()
    except:
        output += [None, None, None]
        if return_inv:
            output += [None, None, None]
    return output


if __name__ == '__main__':
    manga_id = '1-284428'
    filepath = '/Volumes/Work/GZ3D/MPL5_fits'
    filename = '1-284428_127_5679016.fits.gz'

    c = Cube(mangaid=manga_id)
    gz3d = gz3d_fits('{0}/{1}'.format(filepath, filename))
    gz3d.make_all_spaxel_masks(grid_size=c.data['FLUX'].data.shape[1:])

    mean_bar = stack_spectra(c, gz3d.bar_mask_spaxel)
    mean_spiral = stack_spectra(c, gz3d.spiral_mask_spaxel)
    ylabel = r'Flux $[\rm 10^{-17}\,erg\,s^{-1}\,cm^{-2}\,\AA^{-1}]$'
    ax, fig = mean_bar.plot(return_figure=True, label='bar')
    mean_spiral.plot(figure=fig, label='spiral')
    ax.set_title(c.mangaid)
    ax.set_ylabel(ylabel)
    plt.legend()
    plt.show()

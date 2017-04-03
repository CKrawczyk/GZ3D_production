from marvin import config
from marvin.tools.cube import Cube
from marvin.tools.spectrum import Spectrum
from gz3d_fits import gz3d_fits
import matplotlib.pyplot as plt
import numpy as np

# set marvin to local mode
# downloading the cubes gives faster access to the spectra
config.mode = 'local'
config.download = True


def load_spectrum(cube, y, x):
    # a helper function to pull the spectra from the cube
    # this is because the built-in method is very slow
    cube_hdu = cube.data
    kwargs = {
        'units': '1E-17 erg/s/cm^2/Ang/spaxel',
        'wavelength': cube_hdu['WAVE'].data,
        'wavelength_unit': 'Angstrom',
        'ivar': cube_hdu['IVAR'].data[:, y, x],
        'mask': cube_hdu['MASK'].data[:, y, x]
    }
    return Spectrum(cube_hdu['FLUX'].data[:, y, x], **kwargs)


def mean_spectrum(spectra, weights=None):
    if weights is None:
        weights = [1] * len(spectra)
    flux = np.array([s.flux for s in spectra])
    ivar = np.array([s.ivar for s in spectra])
    mask = np.array([s.mask for s in spectra])
    weights = np.array([weights]).T
    stack_flux = np.sum(flux * weights, axis=0) / weights.sum()
    # I think this is the correct ivar propagation
    stack_ivar = (weights.sum()**2) / np.sum((weights)**2 / ivar, axis=0)
    # I assume we want a bitwise_or for the masks
    # can also use bitwise_and if that is needed
    stack_mask = np.bitwise_or.reduce(mask, axis=0)
    kwargs = {
        'units': '1E-17 erg/s/cm^2/Ang',
        'wavelength': spectra[0].wavelength,
        'wavelength_unit': 'Angstrom',
        'ivar': stack_ivar,
        'mask': stack_mask
    }
    return Spectrum(stack_flux, **kwargs)


def stack_spectra(cube, maks):
    yy, xx = np.where(mask > 0)
    if len(yy) > 0:
        mask_value = mask[yy, xx]
        spectra = [load_spectrum(cube, y, x) for y, x in zip(yy, xx)]
        return mean_spectrum(spectra, weights=mask_value)
    else:
        return None


if __name__ == '__main__':
    manga_id = '1-284428'
    filepath = '/Volumes/Work/GZ3D/MPL5_fits'
    filename = '1-284428_127_5679016.fits.gz'

    c = Cube(mangaid=manga_id)
    gz3d = gz3d_fits('{0}/{1}'.format(filepath, filename))
    gz3d.make_all_spaxel_masks()

    mean_bar = stack_spectra(c, gz3d.bar_mask_spaxel)
    mean_spiral = stack_spectra(c, gz3d.spiral_mask_spaxel)
    mean_bar.plot()
    mean_spiral.plot()
    plt.show()

from marvin.tools.spectrum import Spectrum
from marvin.tools.cube import Cube
import numpy as np


# make a subclass of Cube to give faster loding of spaxels
# I only need the spectrum, so no properties need to be loaded
class CubeFast(Cube):
    def __getitem__(self, xy):
        """Returns the spaxel for ``(x, y)`` without propserties"""
        return self.getSpaxel(x=xy[0], y=xy[1], properties=False, xyorig='lower')


# subclass Spectrum to add basic opperations
# such as +, -, *, /, that handel ivar and bit masks correctly
class SpectrumStacker(Spectrum):
    def inside_ifu(self):
        return not any((2**0) & self.mask)

    def _checkUnits(self, s):
        if self.units != s.units:
            raise ValueError('Units must match')

    def _checkWavelength(self, s):
        if all(self.wavelength != s.wavelength):
            raise ValueError('Wavelength must match')

    def _ivar_sum(self, s):
        var1 = 1.0 / self.ivar
        var2 = 1.0 / s.ivar
        return 1 / (var1 + var2)

    def __add__(self, s):
        self._checkUnits(s)
        self._checkWavelength(s)
        flux = self.flux + s.flux
        ivar = self._ivar_sum(s)
        mask = np.bitwise_or(self.mask, s.mask)
        kwargs = {
            'units': self.units,
            'wavelength': self.wavelength,
            'wavelength_unit': self.wavelength_unit,
            'ivar': ivar,
            'mask': mask
        }
        return SpectrumStacker(flux, **kwargs)

    def __radd__(self, s):
        if s == 0:
            return self
        else:
            return self.__add__(s)

    def __sub__(self, s):
        self._checkUnits(s)
        self._checkWavelength(s)
        flux = self.flux - s.flux
        ivar = self._ivar_sum(s)
        mask = np.bitwise_or(self.mask, s.mask)
        kwargs = {
            'units': self.units,
            'wavelength': self.wavelength,
            'wavelength_unit': self.wavelength_unit,
            'ivar': ivar,
            'mask': mask
        }
        return SpectrumStacker(flux, **kwargs)

    def __mul__(self, n):
        flux = self.flux * n
        ivar = self.ivar / (n**2)
        kwargs = {
            'units': self.units,
            'wavelength': self.wavelength,
            'wavelength_unit': self.wavelength_unit,
            'ivar': ivar,
            'mask': self.mask
        }
        return SpectrumStacker(flux, **kwargs)

    def __rmul__(self, s):
        return self.__mul__(s)

    def __div__(self, n):
        flux = self.flux / np.float(n)
        ivar = self.ivar * (n**2)
        kwargs = {
            'units': self.units,
            'wavelength': self.wavelength,
            'wavelength_unit': self.wavelength_unit,
            'ivar': ivar,
            'mask': self.mask
        }
        return SpectrumStacker(flux, **kwargs)

    def __truediv__(self, n):
        return self.__div__(n)

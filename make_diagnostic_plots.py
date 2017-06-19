from diagnostic_plots import diagnostic_plots
from astropy.table import Table
import os
import progressbar as pb


filepath = '/Volumes/Work/GZ3D/MPL5_fits'
filenames = Table.read('bar_or_spiral.txt', format='ascii.csv')
output_folder = '/Volumes/Work/GZ3D/MPL5_plots_theta'


widgets = ['Plot: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
pbar = pb.ProgressBar(widgets=widgets, max_val=len(filenames))
pbar.start()

for fdx, filename in enumerate(filenames['filename']):
    output_name = os.path.join(output_folder, filename.split('.')[0])
    input_name = os.path.join(filepath, filename)
    if os.path.isfile(output_name):
        pbar.update(fdx + 1)
        continue
    diagnostic_plots(input_name, output_name=output_name, fdx=fdx, oi_sf=False)
    pbar.update(fdx + 1)
pbar.finish()

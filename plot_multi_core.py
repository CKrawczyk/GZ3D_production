from diagnostic_plots import diagnostic_plots
from astropy.table import Table
import os
import progressbar as pb
from multiprocessing.pool import Pool
import multiprocessing as mp


class Counter(object):
    def __init__(self, initval=0):
        self.val = mp.Value('i', initval)
        self.lock = mp.Lock()

    def change_value(self, v):
        with self.lock:
            self.val.value = v

    def increment(self):
        with self.lock:
            self.val.value += 1

    def value(self):
        with self.lock:
            return self.val.value


filepath = '/Volumes/Work/GZ3D/MPL5_fits'
filenames = Table.read('bar_or_spiral.txt', format='ascii.csv')
in_names = [os.path.join(filepath, f) for f in filenames['filename']]

output_folder = '/Volumes/Work/GZ3D/MPL5_plots_theta'
out_names = [os.path.join(output_folder, f.split('.')[0]) for f in filenames['filename']]

fdx = range(len(filenames))
cdx = Counter()

widgets = ['Plot: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
pbar = pb.ProgressBar(widgets=widgets, max_val=len(filenames))
pbar.start()


def wrapper(args):
    filename, output_name, fdx = args
    if os.path.isfile(output_name):
        cdx.increment()
        pbar.update(cdx.value)
        return 'already done'
    diagnostic_plots(filename, output_name=output_name, fdx=fdx)
    cdx.increment()
    pbar.update(cdx.value)
    return 'done'


pool = Pool(processes=3)
pool.map(wrapper, zip(in_names, out_names, fdx))
pool.close()
pbar.finish()

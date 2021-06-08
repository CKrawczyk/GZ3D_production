import glob
import os

from astropy.io import fits
from astropy.table import Table
from tqdm import tqdm

input_folder = '/Volumes/Work/GZ3D/data_run_final_fits_no_user'

output_max_mask_count = os.path.join(input_folder, 'gz3d_max_mask_count.csv')
metadata_counter = 0
fits_list = glob.glob(input_folder + '/*fits.gz')

for f in tqdm(fits_list):
    file_name = os.path.basename(f)[:-3]
    subject_id = file_name.split('_')[-1].split('.')[0]
    hdu_list = fits.open(f)

    row = Table({
        'file_name': [file_name],
        'zooniverse_subject_id': [subject_id],
        'spiral_mask_max_count': [hdu_list[3].data.max()],
        'bar_mask_max_count': [hdu_list[4].data.max()]
    })

    row = row.to_pandas()

    if metadata_counter == 0:
        # make new file
        row.to_csv(output_max_mask_count, mode='w', index=False)
        metadata_counter += 1
    else:
        # append to existing file
        row.to_csv(output_max_mask_count, mode='a', index=False, header=False)

    hdu_list.close()

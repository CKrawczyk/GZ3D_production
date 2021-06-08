import glob
import os

from astropy.io import fits
from astropy.table import Table, vstack, hstack
from subprocess import call
from tqdm import tqdm

input_folder = '/Volumes/Work/GZ3D/Data_run_final_fits'

output_folder = '/Volumes/Work/GZ3D/Data_run_final_fits_no_user'

observed_data_hdu_list = fits.open('/Volumes/Work/GZ3D/FinalMaNGAgalaxies_bundlesize_NSA_HIflag.fits')
observed_data = Table(observed_data_hdu_list[1].data)
# strip whitespace from ids
observed_ids = [id.strip() for id in observed_data['MANGAID_1']]
observed_data_hdu_list.close()

metadata_counter = 0
centers_counter = 0
stars_counter = 0

output_metadata = os.path.join(output_folder, 'gz3d_metadata.csv')
output_centers = os.path.join(output_folder, 'galaxy_centers.csv')
output_stars = os.path.join(output_folder, 'foreground_stars.csv')

fits_list = glob.glob(input_folder + '/*fits.gz')

for f in tqdm(fits_list):
    file_name = os.path.basename(f)[:-3]
    subject_id = file_name.split('_')[-1].split('.')[0]
    hdu_list = fits.open(f)

    # create tables for summary files
    metadata = Table(hdu_list[5].data)
    manga_id = metadata['MANGAID'][0].strip()

    centers = Table(hdu_list[6].data)
    if len(centers) > 0:
        centers.add_column(file_name, name='file_name', index=0)
        centers.remove_columns(['x', 'y', 'var_x', 'var_y', 'var_x_y', 'var_ra', 'var_dec', 'var_ra_dec'])
        centers = centers.to_pandas()
        if centers_counter == 0:
            # make new file
            centers.to_csv(output_centers, mode='w', index=False)
            centers_counter += 1
        else:
            # append to existing file
            centers.to_csv(output_centers, mode='a', index=False, header=False)

    stars = Table(hdu_list[7].data)
    if len(stars) > 0:
        stars.add_column(file_name, name='file_name', index=0)
        stars.remove_columns(['x', 'y', 'var_x', 'var_y', 'var_x_y', 'var_ra', 'var_dec', 'var_ra_dec'])
        stars = stars.to_pandas()
        if stars_counter == 0:
            # make new file
            stars.to_csv(output_stars, mode='w', index=False)
            stars_counter += 1
        else:
            # append to existing file
            stars.to_csv(output_stars, mode='a', index=False, header=False)

    new_metadata_columns = Table({
        'file_name': [file_name],
        'zooniverse_subject_id': [subject_id],
        'number_galaxy_centers': [len(centers)],
        'number_foreground_stars': [len(stars)],
        'observed': [manga_id in observed_ids]
    })

    new_metadata = hstack([new_metadata_columns, metadata])
    new_metadata = new_metadata.to_pandas()

    if metadata_counter == 0:
        # make new file
        new_metadata.to_csv(output_metadata, mode='w', index=False)
        metadata_counter += 1
    else:
        # append to existing file
        new_metadata.to_csv(output_metadata, mode='a', index=False, header=False)

    # strip user_name and user_id
    for hdu_index in [8, 9, 10]:
        raw_table = Table(hdu_list[hdu_index].data)
        raw_table.remove_columns(['user_name', 'user_id'])
        hdu_list[hdu_index] = fits.table_to_hdu(raw_table)

    output_name = os.path.join(output_folder, file_name)
    hdu_list.writeto(output_name)
    # compress the fits file
    call(['gzip', output_name])
    hdu_list.close()

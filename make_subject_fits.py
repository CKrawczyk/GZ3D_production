# this is an all-in-one subject aggregating code for GZ: 3D
# Basic idea:
#   1. Loop over all subjects
#    2. Grab all classifications for that subject (from any workflow)
#    3. Aggregate based on workflow
#    4. Save all information for the subject to a FITS file
#
# FITS file structrue:
#    HDU 0: SDSS cutout image (that the used saw)
#    HDU 1: [image] Pixel mask of clustering resutls for galaxy center(s) ("2 sigma" ellipse about clustered points with value equal to number of points used to make the cluster)
#    HDU 2: [image] Pixel mask of clustering resutls for star(s) ("2 sigma" ellipse about clustered points with value equal to number of points used to make the cluster)
#    HDU 3: [image] Pixel mask of spiral arm location(s) (cleaned classification counts for each pixel)
#    HDU 4: [image] Pixel mask of bar location (cleaned classification counts for each pixel)
#    HDU 5: [table] Image metadata
#    HDU 6: [table] Center cluster data table (in pix and RA-DEC coords)
#    HDU 7: [table] Star cluster data table (in pix and RA-DEC coords)
#    HDU 8: [table] raw center and star classifications
#    HDU 9: [table] raw spiral arm classifications
#    HDU 10: [table] raw bar classifications

from panoptes_client import SubjectSet, Panoptes, Classification
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import progressbar as pb
from collections import OrderedDict
from sklearn.cluster import DBSCAN
from cluster_points import cov_to_ellipse
import numpy as np
import scipy.linalg as sl
import matplotlib.pyplot as plt
from matplotlib.path import Path
from shapely.geometry import LineString
import json
import os.path
from subprocess import call
Panoptes.connect()
widgets = ['Aggregate: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]

# metadata on each subject is in this format
metadata_dtype = [
    ('ra', '>f4'),
    ('dec', '>f4'),
    ('MANGAID', 'S11'),
    ('IAUNAME', 'S19'),
    ('IFUDESIGNSIZE', '>f8'),
    ('#MANGA_TILEID', '>f8'),
    ('NSAID', '>i8'),
    ('explorer_link', 'S90'),
    ('GZ2_total_classifications', '>i2'),
    ('GZ2_bar_votes', '>i2'),
    ('GZ2_spiral_votes', '>i2'),
    ('specobjid', '>i8'),
    ('dr7objid', '>i8'),
    ('dr8objid', '>i8'),
    ('gz2_sample', 'S70')
]


def define_wcs(ra, dec, scale=0.099, size_pix=np.array([525, 525])):
    """
    Given what we know about the scale of the image,
    define a nearly-correct world coordinate system to use with it.
    """
    w = WCS(naxis=2)
    w.wcs.crpix = size_pix / 2
    w.wcs.crval = np.array([ra, dec])
    w.wcs.cd = np.array([[-1, 0], [0, -1]]) * scale / 3600.
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.cunit = ['deg', 'deg']
    return w


def cov_to_ellipse_params(cov, pos, nstd=1):
    eigvec, eigval, V = sl.svd(cov, full_matrices=False)
    theta = np.degrees(np.arctan2(eigvec[1, 0], eigvec[0, 0]))
    a, b = nstd * np.sqrt(eigval)
    return a, b, np.deg2rad(theta)


def inside_ellipse(coords, center, a, b, theta):
    # subtract center
    t_coords = np.array(coords) - center
    # rotate coords by theta
    R = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
    t_coords = np.dot(t_coords, R)
    # divide by semi- axies
    t_coords /= np.array([a, b])
    # in this coord system the ellipse is a unit circle
    # check where x**2 + y**2 <= 1
    inside = (t_coords**2).sum(axis=1) <= 1
    return inside


def make_cluster_table():
    return OrderedDict([
        ('x', []),
        ('y', []),
        ('var_x', []),
        ('var_y', []),
        ('var_x_y', []),
        ('ra', []),
        ('dec', []),
        ('var_ra', []),
        ('var_dec', []),
        ('var_ra_dec', []),
        ('count', [])
    ])


def cluster(X, dimensions, coords, wcs):
    mask = np.zeros(dimensions)
    db = DBSCAN(eps=5, min_samples=3).fit(X)
    cluster_table = make_cluster_table()
    X_ra_dec = wcs.wcs_pix2world(X, 1)
    for k in set(db.labels_):
        if k > -1:
            idx = db.labels_ == k
            cluster_table['count'].append(idx.sum())
            m = X[idx].mean(axis=0)
            cluster_table['x'].append(m[0])
            cluster_table['y'].append(m[1])
            cov = np.cov(X[idx].T)
            cluster_table['var_x'].append(cov[0, 0])
            cluster_table['var_y'].append(cov[1, 1])
            cluster_table['var_x_y'].append(cov[0, 1])
            a, b, theta = cov_to_ellipse_params(cov, m, nstd=2)
            if (a > 0) and (b > 0):
                inside = inside_ellipse(coords, m, a, b, theta).reshape(*dimensions)
                mask += inside * idx.sum()
            # now calculate center and cov in RA DEC
            m_ra_dec = X_ra_dec[idx].mean(axis=0)
            cluster_table['ra'].append(m_ra_dec[0])
            cluster_table['dec'].append(m_ra_dec[1])
            cov_ra_dec = np.cov(X_ra_dec[idx].T)
            cluster_table['var_ra'].append(cov_ra_dec[0, 0])
            cluster_table['var_dec'].append(cov_ra_dec[1, 1])
            cluster_table['var_ra_dec'].append(cov_ra_dec[0, 1])
    return mask, fits.table_to_hdu(Table(cluster_table))


def path(X, dimensions, coords):
    # X is list of lists
    mask = np.zeros(dimensions)
    # remove self intersecting paths
    for p in X:
        codes = [Path.MOVETO] + [Path.LINETO] * len(p)
        if len(p) > 1:
            if not LineString(p).is_simple:
                continue
        p_closed = p + [p[0]]
        mpl_path = Path(p_closed, codes=codes)
        inside = mpl_path.contains_points(coords).reshape(*dimensions)
        mask += inside
    return mask


def make_classification_table(*args):
    table = OrderedDict([
        ('classification_id', []),
        ('user_id', []),
        ('time_stamp', [])
    ])
    for a in args:
        table[a] = []
    return table


def record_base_classification(c, table):
    table['classification_id'].append(c.id)
    if 'user' in c.raw['links']:
        table['user_id'].append(c.raw['links']['user'])
    else:
        table['user_id'].append('')
    table['time_stamp'].append(c.raw['created_at'])


def mask_process(workflow_id, subject_id, column_name, wcs_header, dimensions, coords):
    classifications_table = make_classification_table(column_name)
    all_paths = []
    for c in Classification.where(scope='project', workflow_id=workflow_id, subject_id=subject_id):
        record_base_classification(c, classifications_table)
        paths = []
        annotations = c.raw['annotations'][0]['value']
        for a in annotations:
            if 'points' in a:
                points = [[p['x'], p['y']] for p in a['points']]
                paths.append(points)
                all_paths.append(points)
        classifications_table[column_name].append(json.dumps(paths))
    table_hdu = fits.table_to_hdu(Table(classifications_table))
    if len(all_paths):
        path_mask = path(all_paths, dimensions, coords)
    else:
        path_mask = np.zeros(dimensions)
    mask_hdu = fits.ImageHDU(data=path_mask, header=wcs_header)
    return table_hdu, mask_hdu


def make_subject_fits(full_subject_set_id, center_workflow_id, spiral_workflow_id, bar_workflow_id, dimensions=[525, 525], image_location='/Volumes/SD_Extra/manga_images_production/MPL5', output='MPL5_fits'):
    blank_mask = np.zeros(dimensions)
    coords = [[x, y] for y in xrange(dimensions[1]) for x in xrange(dimensions[0])]
    fullsample = SubjectSet.find(full_subject_set_id)
    subjects = fullsample.subjects()
    pbar = pb.ProgressBar(widgets=widgets, maxval=subjects.meta['count'])
    pbar.start()
    idx = 0
    for subject in subjects:
        subject_metadata = Table(dtype=metadata_dtype)
        subject_metadata.add_row(tuple(subject.raw['metadata'][key] for key in subject_metadata.dtype.names))
        subject_metadata.rename_column('#MANGA_TILEID', 'MANGA_TILEID')
        subject_metadata_hdu = fits.table_to_hdu(subject_metadata)
        loc = '{0}/{1}_{2}.jpg'.format(image_location, subject_metadata['MANGAID'][0], int(subject_metadata['IFUDESIGNSIZE'][0]))
        output_name = '{0}/{1}_{2}_{3}.fits'.format(output, subject_metadata['MANGAID'][0], int(subject_metadata['IFUDESIGNSIZE'][0]), subject.id)
        if os.path.isfile('{0}.gz'.format(output_name)):
            # don't process the file if it already exists
            idx += 1
            pbar.update(idx)
            continue
        image = plt.imread(loc, format='jpeg')
        wcs = define_wcs(subject_metadata['ra'][0], subject_metadata['dec'][0])
        wcs_header = wcs.to_header()
        orig_image_hdu = fits.PrimaryHDU(data=image, header=wcs_header)
        # process data from center(s) and star(s) points
        center_classifications = make_classification_table('center_points', 'star_points')
        all_center = []
        all_star = []
        for c in Classification.where(scope='project', workflow_id=center_workflow_id, subject_id=subject.id):
            record_base_classification(c, center_classifications)
            center_points = []
            star_points = []
            points = c.raw['annotations'][0]['value']
            for p in points:
                if ('x' in p) and ('y' in p):
                    if p['tool'] == 0:
                        # somehow the workflow_id got messed up for some classifications
                        # so a try statement is needed
                        loc = [p['x'], p['y']]
                        center_points.append(loc)
                        all_center.append(loc)
                    elif p['tool'] == 1:
                        loc = [p['x'], p['y']]
                        star_points.append(loc)
                        all_star.append(loc)
            center_classifications['center_points'].append(json.dumps(center_points))
            center_classifications['star_points'].append(json.dumps(star_points))
        center_star_table_hdu = fits.table_to_hdu(Table(center_classifications))
        # cluster points and make image masks
        if len(all_center):
            center_mask, center_table_hdu = cluster(np.array(all_center), dimensions, coords, wcs)
        else:
            center_mask = blank_mask
            center_table_hdu = fits.table_to_hdu(Table(make_cluster_table()))
        center_hdu = fits.ImageHDU(data=center_mask, header=wcs_header)
        if len(all_star):
            star_mask, star_table_hdu = cluster(np.array(all_star), dimensions, coords, wcs)
        else:
            star_mask = blank_mask
            star_table_hdu = fits.table_to_hdu(Table(make_cluster_table()))
        star_hdu = fits.ImageHDU(data=star_mask, header=wcs_header)
        # spiral arms
        spiral_table_hdu, spiral_hdu = mask_process(spiral_workflow_id, subject.id, 'spiral_paths', wcs_header, dimensions, coords)
        # bars
        bar_table_hdu, bar_hdu = mask_process(bar_workflow_id, subject.id, 'bar_paths', wcs_header, dimensions, coords)
        # make fits file
        hdu_list = fits.HDUList([orig_image_hdu, center_hdu, star_hdu, spiral_hdu, bar_hdu, subject_metadata_hdu, center_table_hdu, star_table_hdu, center_star_table_hdu, spiral_table_hdu, bar_table_hdu])
        hdu_list.writeto(output_name)
        # compress the fits file
        call(['gzip', output_name])
        # update progressbar
        idx += 1
        pbar.update(idx)
    pbar.finish()


if __name__ == '__main__':
    make_subject_fits(8199, 3513, 3515, 3514, dimensions=[525, 525], image_location='/Volumes/SD_Extra/manga_images_production/MPL5', output='/Volumes/Work/GZ3D/MPL5_fits')

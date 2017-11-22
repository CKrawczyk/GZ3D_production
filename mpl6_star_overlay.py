import pandas
from panoptes_aggregation.csv_utils import unjson_dataframe
from astropy.wcs import WCS
import numpy as np
from gz3d_fits import cov_to_ellipse
from plot_fits_files import set_up_axes
from make_subject_fits import define_wcs
import matplotlib.pyplot as plt
import mpl_style
import progressbar as pb

plt.style.use(mpl_style.style1)

path = '/Volumes/Work/GZ3D'

data = pandas.read_csv(path + '/point_reducer_galaxy_and_star_mpl6.csv')
unjson_dataframe(data)

widgets = ['Plot: ', pb.Percentage(), ' ', pb.Bar(marker='0', left='[', right=']'), ' ', pb.ETA()]
pbar = pb.ProgressBar(widgets=widgets, max_value=len(data))

subjects = pandas.read_csv(path + '/galaxy-zoo-3d-subjects-mpl5-mpl6.csv')
subjects.columns = [
    'subject_id',
    'project_id',
    'workflow_id',
    'subject_set_id',
    'data.metadata',
    'data.locations',
    'classifications_count',
    'retired_at',
    'retirement_reason'
]
unjson_dataframe(subjects)

pbar.start()
for index, reduction in data.iterrows():
    if (not np.isnan(reduction['data.T0_tool1_clusters_count']).any()) and (max(reduction['data.T0_tool1_clusters_count']) >= 10):
        for idx, count in enumerate(reduction['data.T0_tool1_clusters_count']):
            x_list = []
            y_list = []
            ellip_list = []
            if count >= 10:
                x = reduction['data.T0_tool1_clusters_x'][idx]
                y = reduction['data.T0_tool1_clusters_y'][idx]
                x_list.append(x)
                y_list.append(y)
                pos = np.array([x, y])
                cov = np.array([
                    [reduction['data.T0_tool1_clusters_var_x'][idx], reduction['data.T0_tool1_clusters_var_x_y'][idx]],
                    [reduction['data.T0_tool1_clusters_var_x_y'][idx], reduction['data.T0_tool1_clusters_var_y'][idx]]
                ])
                ellip_list.append(cov_to_ellipse(cov, pos, nstd=2, edgecolor='C3', facecolor='none', lw=2, alpha=0.75))
        # plot original image
        sdx = (subjects.subject_id == reduction.subject_id) & (subjects.workflow_id == reduction.workflow_id) & (subjects.subject_set_id == 8215)
        subject_metadata = subjects['data.metadata'].loc[sdx].item()
        wcs = define_wcs(float(subject_metadata['ra']), float(subject_metadata['dec']))
        manga_name = '{0}_{1}'.format(subject_metadata['MANGAID'], int(float(subject_metadata['IFUDESIGNSIZE'])))
        loc = '/Volumes/SD_Extra/manga_images_production/non_MPL5/{0}.jpg'.format(manga_name)
        image = plt.imread(loc, format='jpeg')
        fig = plt.figure(1)
        ax = plt.subplot(111, projection=wcs)
        set_up_axes(ax, color_grid='C7')
        ax.imshow(image)
        ax.plot(x, y, 'x', color='C3')
        for e in ellip_list:
            ax.add_artist(e)
        plt.savefig('{0}/MPL6_stars/{1}.png'.format(path, manga_name))
        plt.close(fig)
    pbar.update(index + 1)
pbar.finish()

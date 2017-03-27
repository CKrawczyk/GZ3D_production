## The data aggregating and plotting code of [Galaxy Zoo: 3D](https://www.zooniverse.org/projects/klmasters/galaxy-zoo-3d)

+ `make_subject_fits.py`: This code uses the `panoptes-python-client` to grab and aggregate the classifications for a point drawing workflow, a polygon drawing workflow, and a poly-freehand drawing workflow.  All classification data and resulting pixel masks are stored in a `.fits` file.
  + The output of the point aggregation is any cluster of points with more than three classifications within 5 pixels of each other (using `DBSCAN`).  The center and 2-sigma extent of each cluster is found and a pixel mask is created for these ellipses.
  + The output of the polygon and poly-freehand aggregation is a pixel mask a count of the number of classifications that contained that pixel.  It also ignores any self-intersecting polygons (as they don't have a well defined "inside").
+ `mpl_style.py`: A custom style file for `matplotlib`
+ `plot_fits_files.py`: This code takes the `.fits` files created above and makes plots showing the various pixel masks
  + Example plot: ![1-167242_127_5679242.png](https://panoptes-uploads.zooniverse.org/production/project_attached_image/c0d94926-efa7-4ea3-80a1-a92e86b852b3.png)

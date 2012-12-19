# setup scipy
# setup esutil -r ~esheldon/exports/esutil-local

# check  ~esheldon/svn/esutil/esutil/wcsutil.py to see how the variables are used

import cPickle as pickle
# wcs_example_file='/astro/u/esheldon/oh/desfiles/wcsexample/wcs-decam--24--15-i-6.pickle'
wcs_example_file='wcs-decam--24--15-i-6.pickle'
many_wcs = pickle.load(open(wcs_example_file))

# print many_wcs[1]

for ccd in many_wcs:
    print many_wcs[ccd]

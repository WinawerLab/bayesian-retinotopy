####################################################################################################
# v123/__init__.py
# This file defines all the basic data-loading and analysis functions used in the accompanying
# iPython notebook.
# by Noah C. Benson

# Import the basics:
import numpy as np
import scipy as sp
import nibabel
import os, sys, math, warnings, time

# Warnings are errors!
warnings.filterwarnings('error')

# Now, import Neuropythy, our core library:
import neuropythy as ny

# This is the configuration hash; it tells us where to find data files and should be edited when
# running these analyses to point to the correct directories.
_configuration = {
    # The data_root is where to find the data downloaded from the Open Science Foundation project
    # page associated with this repository. See the README file in the repository root for more
    # information. If you use the download.sh script, this will be <__file__>/../data_root.
    'data_root': None
}

subject_names = {}
subject_datafiles = {}
subject_dataset_count = {}
subject_datasets = {}

def data_root(*args):
    '''
    v123.data_root() yields the current data root directory for the v123 analysis project.
    v123.data_root(root) sets the data root to be the given root directory and resets all cache
      that might now be invalid.
    '''
    global subject_names, subject_datafiles, subject_dataset_count, subject_datasets
    if not args:
        dr = _configuration['data_root']
        if dr is None or not os.path.exists(dr) or not os.path.isdir(dr):
            raise ValueError('Cannot find data_root: %s' % dr)
    else:
        if len(args) > 1: dr = os.path.join(*args)
        else: dr = args[0]
        _configuration['data_root'] = dr
        clear_cache()
        if not (os.path.exists(dr) and os.path.isdir(dr)): return dr
        fspath = freesurfer_path()
        # We want to run extra initialization if this is a real path:
        if dr is not None and os.path.exists(dr) and os.path.isdir(dr):
            # Setup freesurfer; make sure the directory is on the path:
            fspath = freesurfer_path()
            if os.path.exists(fspath) and os.path.isdir(fspath) \
               and fspath not in ny.freesurfer.subject_paths():
                ny.freesurfer.add_subject_path(fspath)
            # Get our subject database put together
            subject_names = [s for s in os.listdir(retinotopy_path()) if s[0] == 'S']
            # Get the datasets for each subject as well;
            # NOTE: 0 is always the gold-standard dataset; the remainder are partial datasets
            subject_datafiles = {
                sub: {
                    hem: {
                        meas: {
                            ds: os.path.join(retinotopy_path(), sub, flnm)
                            for flnm in dirfiles
                            if flnm[0] == hem[0] and flnm[-9:-4] == meas
                            for ds in [int(flnm[3:5])]}
                        for meas in ['angle', 'eccen', 'vexpl', 'prfsz']}
                    for hem in ['lh', 'rh']}
                for sub in subject_names
                for dirfiles in [os.listdir(os.path.join(retinotopy_path(), sub))]}
            subject_dataset_count = {
                sub: len(subject_datafiles[sub]['lh']['angle'])
                for sub in subject_names}
            subject_datasets = {
                sub: subject_datafiles[sub]['lh']['angle'].keys()
                for sub in subject_names}
    return dr
def retinotopy_path():
    '''
    v123.retinotopy_path() yields the current directory containing the v123 analysis project's
      individual subject retinotopy surface-overlay data.
    '''
    return os.path.join(data_root(), 'retinotopy')
def cache_path():
    '''
    v123.cache_path() yields the current directory containing the v123 analysis project's
      individual subject data cache, where intermediate files are stored for faster loading.
    '''
    return os.path.join(data_root(), 'cache')
def freesurfer_path():
    '''
    v123.freesurfer_path() yields the current directory for the v123 analysis project's FreeSurfer
      subjets.
    '''
    return os.path.join(data_root(), 'freesurfer_subjects')
def analyses_path():
    '''
    v123.analyses_path() yields the current directory for the v123 analysis project's results;
      outputs of analyses should be written here.
    '''
    return os.path.join(data_root(), 'analyses')

# Create/load the Model of V1/V2/V3 that we will be using
v123_model = ny.V123_model()

# This function converts a variance-explained measurement into a weight by
# running it through an error function:
def vexpl_to_weight(vexpl):
    ierf05 = sp.special.erfinv(0.5)
    vexpl = np.asarray(vexpl)
    w = vexpl[vexpl > 0.2]
    return sp.special.erf((ierf05/np.median(w)) * (~np.isclose(vexpl, 0)) * vexpl)

# The following functions are for extracting data from the above datafiles;
# all of these cache the data as they go so that data is not reloaded more
# that once per session.
_subject_data_cache = {}
_measure_names = ['polar_angle', 'eccentricity', 'variance_explained', 'weight']
def subject_data(sub, hem, ds, meas=None):
    global _subject_data_cache
    from nibabel.freesurfer.mghformat import load as mghload
    if isinstance(sub, (int, long)): sub = ('S%04d' % sub)
    if sub not in subject_names:
        raise ValueError('subject %s not found in database' % sub)
    hem = hem.lower()
    if hem != 'lh' and hem != 'rh':
        raise ValueError('hem must be lh or rh')
    if sub == 'aggregate' and hem == 'rh':
        raise ValueError('aggregate subject has only a left-hemisphere')
    if meas is not None:
        meas = meas.lower()
        if meas in ['angle', 'eccen', 'vexpl', 'prfsz']:
            meas = {'angle': 'polar_angle'       ,
                    'eccen': 'eccentricity'      ,
                    'vexpl': 'variance_explained',
                    'prfsz': 'prf_size'          }[meas]
        if meas not in _measure_names:
            raise ValueError('meas type not valid: %s' % meas)
    if ds not in subject_datasets[sub]:
        raise ValueError('subject %s does not have a dataset %d' % (sub, ds))
    if sub not in _subject_data_cache: _subject_data_cache[sub] = {}
    subdat = _subject_data_cache[sub]
    if hem not in subdat: subdat[hem] = {}
    hemdat = subdat[hem]
    if ds not in hemdat:
        # now we load the data...
        fs = subject_datafiles[sub][hem]
        hemdat[ds] = {
            'polar_angle':        mghload(fs['angle'][ds]).get_data().flatten(),
            'eccentricity':       mghload(fs['eccen'][ds]).get_data().flatten(),
            'variance_explained': mghload(fs['vexpl'][ds]).get_data().flatten(),
            'prf_size':           mghload(fs['prfsz'][ds]).get_data().flatten()}
        hemdat[ds]['weight'] = vexpl_to_weight(hemdat[ds]['variance_explained'])
    dsdat = hemdat[ds]
    return dsdat if meas is None else dsdat[meas]
def aggregate_data(meas=None):
    global _subject_data_cache
    from nibabel.freesurfer.mghformat import load as mghload
    if meas is not None:
        meas = meas.lower()
        if meas in ['angle', 'eccen', 'vexpl']:
            meas = {'angle': 'polar_angle'       ,
                    'eccen': 'eccentricity'      ,
                    'vexpl': 'variance_explained',
                    'prfsz': 'prf_size'          }[meas]
        if meas not in _measure_names:
            raise ValueError('meas type not valid: %s' % meas)
    if 'agg' not in _subject_data_cache:
        # now we load the data...
        aggdat = {
            k: mghload(os.path.join(retinotopy_path(), 'aggregate', v)).get_data().flatten()
            for (k,v) in {'polar_angle':  'angle.mgz',
                          'eccentricity': 'eccen.mgz',
                          'weight':       'weight.mgz'}.iteritems()}
        _subject_data_cache['agg'] = aggdat
    return _subject_data_cache['agg']

_subject_hemi_cache = {}
_subject_hemi_cache_ds = {}
def subject_hemi(sub, hem, ds=None):
    global _subject_hemi_cache, _subject_hemi_cache_ds
    if ds is None:
        if sub not in _subject_hemi_cache:
            _subject_hemi_cache[sub] = {}
        scache = _subject_hemi_cache[sub]
        hem = hem.lower()
        if hem not in scache:
            gsdat = subject_data(sub, hem, 0)
            hemi = getattr(ny.freesurfer_subject(sub), hem.upper())
            hemi = hemi.using(
                properties = reduce(
                    lambda p,x: p.without(x) if x in p else p,
                    ['PRF_polar_angle', 'PRF_eccentricity', 'PRF_size',
                     'PRF_variance_explained',
                     'predicted_polar_angle', 'predicted_visual_area',
                     'predicted_eccentricity',
                     'polar_angle', 'eccentricity', 'visual_area'],
                    hemi.properties))
            hemi = hemi.using(
                properties=hemi.properties.using(
                    gold_polar_angle        = gsdat['polar_angle'],
                    gold_eccentricity       = gsdat['eccentricity'],
                    gold_variance_explained = gsdat['variance_explained'],
                    gold_weight             = gsdat['weight'],
                    gold_prf_size           = gsdat['prf_size']))
            scache[hem] = hemi
        return scache[hem]
    elif (sub,hem,ds) in _subject_hemi_cache_ds:
        return _subject_hemi_cache_ds[(sub,hem,ds)]
    else:
        shemi = subject_hemi(sub, hem)
        sdat = subject_data(sub, hem, ds)
        shemi = shemi.using(
            properties=shemi.properties.using(
                    polar_angle        = sdat['polar_angle'],
                    eccentricity       = sdat['eccentricity'],
                    variance_explained = sdat['variance_explained'],
                    weight             = sdat['weight'],
                    prf_size           = sdat['prf_size']))
        _subject_hemi_cache_ds[(sub,hem,ds)] = shemi
        return shemi
def aggregate_hemi():
    global _subject_hemi_cache
    if 'agg' in _subject_hemi_cache: return _subject_hemi_cache['agg']
    aggdat = aggregate_data()
    agghem = ny.freesurfer_subject('fsaverage_sym').LH
    agghem = agghem.using(
        properties = reduce(
            lambda p,x: p.without(x) if x in p else p,
            ['PRF_polar_angle', 'PRF_eccentricity', 'PRF_size',
             'PRF_variance_explained',
             'predicted_polar_angle', 'predicted_visual_area',
             'predicted_eccentricity',
             'polar_angle', 'eccentricity', 'visual_area'],
             agghem.properties))
    agghem = agghem.using(
        properties = agghem.properties.using(
            polar_angle=aggdat['polar_angle'],
            eccentricity=aggdat['eccentricity'],
            weight=aggdat['weight']))
    _subject_hemi_cache['agg'] = agghem
    return agghem
    
_subject_prep_cache = {}
def subject_prep(sub, hemi, ds):
    '''
    subject_prep(sub, hemi, dataset_no) yields a map of data as prepared by the
      neuropythy.vision.register_retinotopy_initialize function; this data should
      be ready for registration.
    '''
    global _subject_prep_cache
    tpl = (sub,hemi,ds)
    if tpl in _subject_prep_cache:
        return _subject_prep_cache[tpl]
    # We need to get the weights right
    p = ny.vision.register_retinotopy_initialize(subject_hemi(sub,hemi,ds),
                                                 v123_model,
                                                 weight_cutoff=0.1)
    _subject_prep_cache[tpl] = p
    return p
def aggregate_prep():
    '''
    aggregate_prep() is like subject_prep but yields a preparation for the
    aggregate dataset instead of an individual subject.
    '''
    global _subject_prep_cache
    if 'agg' in _subject_prep_cache: return _subject_prep_cache['agg']
    p = ny.vision.register_retinotopy_initialize(aggregate_hemi(), v123_model,
                                                 prior=None, resample=None,
                                                 weight_cutoff=0.1,
                                                 partial_voluming_correction=False)
    _subject_prep_cache['agg'] = p
    return p

def auto_cache(filename, calc_fn, cache_directory=Ellipsis, create_dirs=True):
    '''
    auto_cache(filename, calc_fn) first checks to see if filename exists, and, if so,
      unpickles it and returns it; otherwise, runs calc_fn(), pickles the result in
      filename, then returns it.

    The following options are accepted:
      * cache_directory (default: v123.cache_path()) specifies the directory in which
        cache files are saved; the filename argument is a relative path from this
        directory. If this directory does not exist, it is created.
      * create_dirs (default: True) specifies whether or not to automatically create
        directories when they do not exist.
    '''
    import pickle
    if cache_directory is Ellipsis: cache_directory = cache_path()
    # make sure cache path exists
    if not os.path.exists(cache_directory) \
       or (not os.path.isdir(cache_directory) and create_dirs):
        os.makedirs(cache_directory)
    if not os.path.isabs(filename):
        filename = os.path.join(cache_directory, filename)
    filename = os.path.abspath(filename)
    if os.path.exists(filename):
        with open(filename, 'rb') as fl:
            res = pickle.load(fl)
    else:
        res = calc_fn()
        fdir = os.path.join('/', *(filename.split(os.sep)[:-1]))
        if len(fdir) > 0 and not os.path.isdir(fdir) and create_dirs:
            os.makedirs(fdir)
        with open(filename, 'wb') as fl:
            pickle.dump(res, fl)
    return res

def _register_calc_fn(prepfn, steps, scale, ethresh):
    def _calc():
        # Prep the aggregate data
        dat = prepfn()
        # Generate a field
        # We don't use the mesh_field currently, but maybe it will be good for other areas.
        #mesh_field = ny.vision.retinotopy.retinotopy_mesh_field(
        #  dat['map'], dat['model'],
        #  scale=scale, exclusion_threshold=ethresh, 
        #  max_eccentricity=1.0, max_polar_angle=18.0, sigma=6.0)
        anchor_field = ny.vision.retinotopy_anchors(dat['map'], dat['model'],
                                                    scale=scale,
                                                    weight_cutoff=0)
        # Register the data
        reg = ny.registration.mesh_register(
            dat['map'],
            [['edge',      'harmonic',      'scale',1.0],
             ['angle',     'infinite-well', 'scale',1.0],
             ['perimeter', 'harmonic'                  ],
             anchor_field],
            method='random',
            max_steps=steps,
            max_step_size=0.05)
        # and postprocess the registration
        postproc = dat['postprocess_function'](reg)
        tmp = {k:postproc[k]
               for k in ['registered_coordinates', 'prediction', 
                         'initial_polar_angle', 'initial_eccentricity', 'initial_weight',
                         'sub_polar_angle', 'sub_eccentricity', 'sub_weight']}
        tmp['resampled_registered_coordinates'] = postproc['registered_coordinates']
        tmp['resampled_registered_coordinates_3D'] = postproc['finished_registration'].coordinates.T
        tmp['registered_coordinates'] = postproc['registration'].coordinates.T
        return tmp
    return _calc

def aggregate_register(steps=2000, scale=1.0, exclusion_threshold=None):
    '''
    aggregate_register() yields a dictionary of data that is the result of
      registering the 2D mesh constructed in aggregate_prep() to the V1/2/3 model.
      The result is cached so that it is not recalculated in the future.

    The following options are accepted:
      * steps (default: 500) specifies the number of steps to run in the minimization.
      * scale (default: 0.1) specifies the scale of the retinotopy force field term.
    '''
    thresh_str = 'none' if exclusion_threshold is None else \
                 ('%06.3f' % exclusion_threshold)
    flnm = 'steps=%05d_scale=%06.3f_thresh=%s.p' % (steps, scale, thresh_str)
    return auto_cache(
        os.path.join('aggregate', flnm),
        _register_calc_fn(aggregate_prep, steps, scale,
                          exclusion_threshold))

_agg_cache = {}
def aggregate(steps=2000, scale=1.0, exclusion_threshold=None):
    '''
    aggregate() yields a left hemisphere mesh object for the fsaverage_sym subject with
      data from both the group average retinotopy and the 'Benson14' registered
      predictions of retinotopy stored in the properties. The mesh's coordinates
      have been registered to the V123 retinotopy model.
    '''
    global _agg_cache
    tpl = (steps, scale, exclusion_threshold)
    if tpl in _agg_cache: return _agg_cache[tpl]
    dat = aggregate_register(steps=steps, scale=scale,
                             exclusion_threshold=exclusion_threshold)
    pre = dat['prediction']
    mesh = ny.freesurfer_subject('fsaverage_sym').LH.sphere_surface
    mesh = mesh.using(
        coordinates=dat['registered_coordinates'],
        properties=mesh.properties.using(
            PRF_polar_angle=dat['sub_polar_angle'],
            PRF_eccentricity=dat['sub_eccentricity'],
            weight=dat['sub_weight'],
            predicted_polar_angle=pre['polar_angle'],
            predicted_eccentricity=pre['eccentricity'],
            predicted_visual_area=pre['V123_label']))
    _agg_cache[tpl] = mesh
    return mesh

def save_aggregate(directory=None, steps=2000, scale=1.0, create_directory=True):
    '''
    save_aggregate() saves the aggregate data in four files placed in the 
      <analyses directory>/aggregate directory; these files are called:
       * lh.retinotopy.steps=<steps>_scale=<scale>.sphere.reg
       * lh.predict_angle.steps=<steps>_scale=<scale>.mgz
       * lh.predict_eccen.steps=<steps>_scale=<scale>.mgz
       * lh.predict_v123roi.steps=<steps>_scale=<scale>.mgz
    The directory into which the files are saved can be set directly with the directory option. The
    options steps and scale are also accepted.
    '''
    agg = aggregate(steps=steps, scale=scale)
    dr = directory if directory is not None else os.path.join(analyses_path(), 'aggregate')
    if not os.path.exists(dr) and create_directory:
        os.makedirs(dr)
    if not os.path.exists(dr) or not os.path.isdir(dr):
        raise ValueError('Could not create directory: %d' % dr)
    # First, save out the predicted surface overlays:
    _surf_mgh = lambda dat, dt: nibabel.freesurfer.mghformat.MGHImage(
        np.asarray([[dat]], dtype=dt),
        np.eye(4))
    flnm_tag = 'steps=%05d_scale=%05.2f' % (steps, scale)
    flnm_pre_tmpl = 'lh.predict_%s.' + flnm_tag + '.mgz'
    img = _surf_mgh(agg.prop('predicted_polar_angle'), np.float32)
    img.to_filename(os.path.join(dr, flnm_pre_tmpl % 'angle'))
    img = _surf_mgh(agg.prop('predicted_eccentricity'), np.float32)
    img.to_filename(os.path.join(dr, flnm_pre_tmpl % 'eccen'))
    img = _surf_mgh(agg.prop('predicted_visual_area'), np.int32)
    img.to_filename(os.path.join(dr, flnm_pre_tmpl % 'v123roi'))
    # Then, save out registration sphere
    flnm_sph = os.path.join(dr, 'lh.retinotopy.%s.sphere.reg' % flnm_tag)
    nibabel.freesurfer.write_geometry(flnm_sph, agg.coordinates.T, agg.indexed_faces.T)
    return None


def subject_register(sub, hem, ds, steps=2000, scale=1.0, exclusion_threshold=None):
    '''
    subject_register(sub, hem, ds) yields a dictionary of data that is the result
      of registering the 2D mesh constructed in subject_prep(sub, hem, ds) to the
      V1/2/3 model. The result is cached so that it is not recalculated in the 
      future.

    The following options are accepted:
      * steps (default: 500) specifies the number of steps to run in the minimization.
      * scale (default: 0.1) specifies the scale of the retinotopy force field term.
    '''
    thresh_str = 'none' if exclusion_threshold is None else \
                 ('%06.3f' % exclusion_threshold)
    flnm = '%s.ds=%02d_steps=%05d_scale=%06.3f_thresh=%s.p' % (
        hem.lower(), ds, 
        steps, scale,
        thresh_str)
    return auto_cache(
        os.path.join(sub, flnm),
        _register_calc_fn(lambda:subject_prep(sub, hem, ds), steps, scale, 
                          exclusion_threshold))

_sub_cache = {}
def subject(sub, hem, ds, steps=2000, scale=1.0, exclusion_threshold=None):
    '''
    subject(sub, hem, ds) yields a appropriate mesh object for the subject whose subject id is
      given (sub) with data from the appropriate dataset (ds) applied as the properties
      'PRF_polar_angle', 'PRF_eccentricity', and 'weight' as well as the predicted retinotopy, as
      deduced via registration, under the property names 'predicted_polar_angle', 
      'predicted_eccentricity', and 'predicted_visual_area'.  The mesh's coordinates have been
      registered to the V123 retinotopy model, so right hemispheres will be returned as an RHX
      mesh object.
    '''
    global _sub_cache
    tpl = (sub, hem, ds, steps, scale, exclusion_threshold)
    if tpl in _sub_cache: return _sub_cache[tpl]
    dat = subject_register(sub, hem, ds, steps=steps, scale=scale,
                           exclusion_threshold=exclusion_threshold)
    pre = dat['prediction']
    hemi = subject_hemi(sub, hem, ds)
    mesh = hemi.sym_sphere_surface.using(
        coordinates=dat['registered_coordinates'],
        properties=hemi.properties.using(
            PRF_polar_angle=dat['sub_polar_angle'],
            PRF_eccentricity=dat['sub_eccentricity'],
            weight=dat['sub_weight'],
            predicted_polar_angle=pre['polar_angle'],
            predicted_eccentricity=pre['eccentricity'],
            predicted_visual_area=pre['V123_label']))
    _sub_cache[tpl] = mesh
    return mesh

def save_subject(sub, hem, ds, directory=None, create_directory=True):
    '''
    save_subject(sub, hem, ds) saves the provided subject's registration data and predictions to the
      <analyses directory>/<subject name> directory; these files are called:
       * lh.retinotopy.<ds>.sphere.reg OR rhx.retinotopy.<ds>.sphere.reg
       * ?h.predict_angle.<ds>.mgz
       * ?h.predict_eccen.<ds>.mgz
       * ?h.predict_v123roi.<ds>.mgz
    The directory into which the files are saved can be set directly with the directory option. The
    options steps and scale are also accepted.
    '''
    dat = subject(sub, hem, ds)
    dr = directory if directory is not None else os.path.join(analyses_path(), sub)
    if not os.path.exists(dr) and create_directory:
        os.makedirs(dr)
    if not os.path.exists(dr) or not os.path.isdir(dr):
        raise ValueError('Could not create directory: %d' % dr)
    # First, save out the predicted surface overlays:
    _surf_mgh = lambda dat, dt: nibabel.freesurfer.mghformat.MGHImage(
        np.asarray([[dat]], dtype=dt),
        np.eye(4))
    hmname = 'lh' if hem.lower() == 'lh' else 'rhx'
    dsname = '%02d' % ds
    flnm_pre_tmpl = hmname + '.predict_%s.' + dsname + '.mgz'
    img = _surf_mgh(dat.prop('predicted_polar_angle'), np.float32)
    img.to_filename(os.path.join(dr, flnm_pre_tmpl % 'angle'))
    img = _surf_mgh(dat.prop('predicted_eccentricity'), np.float32)
    img.to_filename(os.path.join(dr, flnm_pre_tmpl % 'eccen'))
    img = _surf_mgh(dat.prop('predicted_visual_area'), np.int32)
    img.to_filename(os.path.join(dr, flnm_pre_tmpl % 'v123roi'))
    # Then, save out registration sphere
    flnm_sph = os.path.join(dr, '%s.retinotopy.%s.sphere.reg' % (hem.lower(), dsname))
    nibabel.freesurfer.write_geometry(flnm_sph, dat.coordinates.T, dat.indexed_faces.T)
    return None

_sub_cmag_cache = {}
def subject_cmag(sub, hem):
    '''
    subject_cmag(sub, hem) calculates and yields the cortical magnification for the given subject's
      given hemisphere. 
    '''
    tpl = (sub, hem.lower())
    if tpl in _sub_cmag_cache: return _sub_cmag_cache[tpl]
    hemi = subject_hemi(sub, hem, 0)
    s = subject(sub, hem, 0)
    vlab = s.prop('predicted_visual_area')
    eccs = s.prop('predicted_eccentricity')
    ang0 = s.prop('predicted_polar_angle')
    angs = np.pi/180.0 * (90 - s.prop('polar_angle'))
    (x,y) = (eccs*np.cos(angs), eccs*np.sin(angs))
    cmag_nei = {}
    cmag_pth = {}
    # We look at both white and pial surfaces:
    for surf in ['white', 'pial']:
        msh = getattr(hemi, surf + '_surface')
        # First, do all the paths; we do these once per area
        for area in [1, 2, 3]:
            # start with eccentricity
            for ecc in [0.31, 0.64, 1.25, 2.5, 5.0, 10.0, 20.0]:
                path = np.asarray(list(reversed(range(-90,90,2))), dtype=np.float) * np.pi/180
                path = ecc * np.asarray((np.cos(path), np.sin(path)))
                cm = ny.vision.path_cortical_magnification(msh, path,
                                                           polar_angle=ang0,
                                                           eccentricity=eccs,
                                                           mask=(vlab== area),
                                                           return_all=True)
                cmag_pth[surf + ('_V%d_ecc=%05.2f' % (area, ecc))]  = [
                    (np.asarray(spth), np.asarray(vpth))
                    for (spth,vpth) in zip(cm[0], cm[1])]
            # then polar angle
            for angd in [0.0, 10.0, 45.0, 90.0, 135.0, 170.0, 180.0]:
                angr = np.pi/180*(90 - angd)
                rmtx = [(np.cos(angr), -np.sin(angr)), (np.sin(angr), np.cos(angr))]
                path = np.asarray(range(0,51,1), dtype=np.float) / 50.0 * 20.0
                path = np.dot(rmtx, np.asarray([path,np.zeros(path.shape)]))
                use_mask = []
                use_angs = []
                # we have to do some odd things to the polar angle/eccen depending on
                # the visual area and angle:
                if (angd == 10.0 or angd == 170.0) and area != 3:
                    # for now, we only do the 'close-to-outer' angles for V3
                    continue
                if angd == 90.0 and area != 1:
                    # we do the 90 degree dorsal or ventral line:
                    nm = surf + ('_%sL_ang=090' % ('D' if area == 2 else 'V'))
                    use_mask = ((vlab == 2) | (vlab == 3))
                    use_mask = use_mask * ((ang0 >= 90.0) if area == 2 else (ang0 <= 90.0))
                    idcs = np.where(vlab == 2)[0]
                    use_angs = np.array(ang0)
                    use_angs[idcs] = 180.0 - use_angs[idcs]
                elif angd == 0.0 or angd == 180.0:
                    # we don't do the outer edges of V3 for now; V2 is done via V1
                    if area == 3 or area == 2: continue
                    nm = surf + ('_V1_ang=%03d' % int(angd))
                    use_mask = ((vlab == 1) | (vlab == 2))
                    use_angs = np.array(ang0)
                    idcs = np.where(vlab == 2)[0]
                    use_angs[idcs] = -use_angs[idcs]
                    use_mask = use_mask * ((angd <= 90.0) if angd == 0.0 else (angd >= 90.0))
                else:
                    nm = surf + ('_V%d_ang=%03d' % (area, int(angd)))
                    use_mask = (vlab == area)
                    use_angs = ang0
                cm = ny.vision.path_cortical_magnification(msh, path,
                                                           polar_angle=use_angs,
                                                           eccentricity=eccs,
                                                           mask=use_mask,
                                                           return_all=True)
                cmag_pth[nm]  = [(np.asarray(spth), np.asarray(vpth))
                                 for (spth,vpth) in zip(cm[0], cm[1])]
        # Then, calculate the neighborhood-based magnification
        cm = ny.vision.neighborhood_cortical_magnification(msh, [x, y])
        for (i,nm) in enumerate(['radial', 'tangential', 'areal']):
            cmag_nei[surf + '_' + nm] = cm[:,i]

    cmag = {'neighborhood': cmag_nei, 'path': cmag_pth}
    _sub_cmag_cache[tpl] = cmag
    return cmag

def save_subject_cmag(sub, hem, directory=None, create_directory=True):
    '''
    save_subject_cmag(sub, hem) saves the data structures found in subject_cmag(sub,hem) out to disk
      in the analyses_path()/<subject name> directory (may be modified with the directory argument).
    '''
    hem = hem.lower()
    s = subject_cmag(sub, hem)
    dr = directory if directory is not None else os.path.join(analyses_path(), sub)
    if not os.path.exists(dr) and create_directory:
        os.makedirs(dr)
    if not os.path.exists(dr) or not os.path.isdir(dr):
        raise ValueError('Could not create directory: %d' % dr)
    # First, save out the predicted surface overlays from the neighborhood datae:
    _surf_mgh = lambda dat, dt: nibabel.freesurfer.mghformat.MGHImage(
        np.asarray([[dat]], dtype=dt),
        np.eye(4))
    flnm_tmpl = hem + '.cm_%s.mgz'
    for (nm,vals) in s['neighborhood'].iteritems():
        _surf_mgh(vals, np.float32).to_filename(os.path.join(dr, flnm_tmpl % nm))
    # Then, save out the paths; this is done in text files
    for (k,cm) in s['path'].iteritems():
        fname = os.path.join(dr, '%s.cmpath_%s.dat' % (hem, k))
        np.savetxt(fname,
                   [(i,a,b,c,d,e)
                    for (i,(spth,vpth)) in zip(range(len(cm)), cm)
                    for ((a,b,c),(d,e)) in zip(spth,vpth)],
                   fmt=('%4d','%8.4f','%8.4f','%8.4f','%8.4f','%8.4f'))
    return None
    
def clear_cache():
    '''
    clear_cache() clears all in-memory caches of the subject database and yields None.
    '''
    global _subject_hemi_cache_ds, _subject_hemi_cache, _subject_data_cache, _subject_prep_cache, \
           _agg_cache
    _subject_hemi_cache_ds = {}
    _subject_hemi_cache = {}
    _subject_data_cache = {}
    _subject_prep_cache = {}
    _subject_cmag_cache = {}
    _agg_cache = {}
    _sub_cache = {}
    return None

# We want to try to setup the data root if possible automatically:
def restore_default_config():
    '''
    v123.restore_default_config() clears the configuration settings of the v123 package and
      resets the cache. If the directory that is selected by default for the data root is found,
      then it is returned, otherwise None is returned.

    By default, the value of the environment variable V123_DATA_ROOT will be used as the root,
    followed by dirname(__file__)/../data_root() (where dirname(__file__) is the v123/ source code
    directory.
    '''
    dr = os.getenv('V123_DATA_ROOT')
    if dr is None or not os.path.exists(dr) or not os.path.isdir(dr):
        dr = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'data_root')
    data_root(dr)
    if dr is not None and os.path.exists(dr) and os.path.isdir(dr):
        return dr
    else:
        return None

# We run this at load-time!
restore_default_config()

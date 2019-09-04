"""Script to run a global SDR analysis."""
import multiprocessing
import math
import glob
import urllib.request
import shutil
import re
import hashlib
import datetime
import os
import zipfile
import logging
import sys

import natcap.invest.sdr
from osgeo import ogr
from osgeo import gdal
from osgeo import osr
import pygeoprocessing
import taskgraph

# set a 1GB limit for the cache
gdal.SetCacheMax(2**30)

WORKSPACE_DIR = 'workspace_global_sdr_dir'
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')
ECOSHARD_DIR = os.path.join(CHURN_DIR, 'ecoshards')
SDR_WORKSPACES_DIR = os.path.join(WORKSPACE_DIR, 'sdr_workspaces')
MOSAIC_DEGREE_CELL_SIZE = 300.0 / 110570
_WGS84_SRS = osr.SpatialReference()
_WGS84_SRS.ImportFromEPSG(4326)
WSGS84_WKT = _WGS84_SRS.ExportToWkt()
DEM_TARGET_NODATA = -32768

N_CPUS = multiprocessing.cpu_count()
TASKGRAPH_REPORTING_FREQUENCY = 5.0
LOGGING_LEVEL = logging.DEBUG

logging.basicConfig(
    level=LOGGING_LEVEL,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

EROSIVITY_URL = r'https://storage.googleapis.com/ecoshard-root/GlobalR_NoPol_compressed_md5_ab6d34ca8827daa3fda42a96b6190ecc.tif'
ERODIBILITY_URL = r'https://storage.googleapis.com/ecoshard-root/pasquale/Kfac_SoilGrid1km_GloSEM_v1.1_md5_e1c74b67ad7fdaf6f69f1f722a5c7dfb.tif'
LULC_URL = r'https://storage.googleapis.com/ipbes-ndr-ecoshard-data/ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_md5_1254d25f937e6d9bdee5779d377c5aa4.tif'
DEM_URL = r'https://storage.googleapis.com/global-invest-sdr-data/global_dem_3s_md5_22d0c3809af491fa09d03002bdf09748.zip'
WATERSHEDS_URL = r'https://storage.googleapis.com/global-invest-sdr-data/watersheds_globe_HydroSHEDS_15arcseconds_md5_c6acf2762123bbd5de605358e733a304.zip'
BIOPHYSICAL_TABLE_URL = r'https://storage.googleapis.com/global-invest-sdr-data/Biophysical_table_ESA_ARIES_RS_md5_e16587ebe01db21034ef94171c76c463.csv'


def main():
    """Entry point."""
    for dir_path in [
            WORKSPACE_DIR, CHURN_DIR, ECOSHARD_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), N_CPUS,
        TASKGRAPH_REPORTING_FREQUENCY)

    root_logger = logging.getLogger()
    root_logger.setLevel(LOGGING_LEVEL)

    lulc_path = os.path.join(ECOSHARD_DIR, os.path.basename(LULC_URL))
    fetch_lulc_task = task_graph.add_task(
        func=url_fetch_and_validate,
        args=(LULC_URL, lulc_path),
        target_path_list=[lulc_path],
        task_name='fetch lulc raster')

    erosivity_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(EROSIVITY_URL))
    fetch_erosivity_task = task_graph.add_task(
        func=url_fetch_and_validate,
        args=(EROSIVITY_URL, erosivity_path),
        target_path_list=[erosivity_path],
        task_name='fetch erosivity raster')

    erodibility_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(ERODIBILITY_URL))
    fetch_erodibility_task = task_graph.add_task(
        func=url_fetch_and_validate,
        args=(ERODIBILITY_URL, erodibility_path),
        target_path_list=[erodibility_path],
        task_name='fetch erodibility raster')

    biophysical_table_path = os.path.join(
        ECOSHARD_DIR, os.path.basename(BIOPHYSICAL_TABLE_URL))
    fetch_biophysical_table_task = task_graph.add_task(
        func=url_fetch_and_validate,
        args=(BIOPHYSICAL_TABLE_URL, biophysical_table_path),
        target_path_list=[biophysical_table_path],
        task_name='fetch biophysical_table raster')

    dem_token_path = os.path.join(
        ECOSHARD_DIR, '%s.COMPLETE' % os.path.basename(DEM_URL))
    fetch_dem_task = task_graph.add_task(
        func=download_validate_and_unzip,
        args=(DEM_URL, ECOSHARD_DIR, dem_token_path),
        target_path_list=[dem_token_path],
        task_name='fetch dem raster')

    watersheds_token_path = os.path.join(
        ECOSHARD_DIR, '%s.COMPLETE' % os.path.basename(WATERSHEDS_URL))
    fetch_watersheds_task = task_graph.add_task(
        func=download_validate_and_unzip,
        args=(WATERSHEDS_URL, ECOSHARD_DIR, watersheds_token_path),
        target_path_list=[watersheds_token_path],
        task_name='fetch watersheds shapefile')

    dem_vrt_path = os.path.join(CHURN_DIR, 'global_dem.vrt')
    dem_vrt_token_path = os.path.join(
        CHURN_DIR, '%s.COMPLETE' % os.path.basename(dem_vrt_path))
    base_raster_pattern = os.path.join(ECOSHARD_DIR, 'global_dem_3s', '*.tif')
    make_dem_task = task_graph.add_task(
        func=make_vrt,
        args=(
            base_raster_pattern, DEM_TARGET_NODATA, dem_vrt_path,
            dem_vrt_token_path),
        dependent_task_list=[fetch_dem_task],
        ignore_path_list=[dem_vrt_path],
        target_path_list=[dem_vrt_token_path],
        task_name='make dem vrt')
    scheduled_watershed_prefixes = set()
    task_graph.join()
    fetch_watersheds_task.join()
    LOGGER.debug('iterating over hydrosheds')
    for watershed_path in glob.glob(os.path.join(
            ECOSHARD_DIR, 'watersheds_globe_HydroSHEDS_15arcseconds',
            '*.shp')):
        LOGGER.debug(watershed_path)
        watershed_basename = os.path.splitext(
            os.path.basename(watershed_path))[0]
        watershed_vector = gdal.OpenEx(watershed_path, gdal.OF_VECTOR)
        watershed_layer = watershed_vector.GetLayer()
        for watershed_feature in watershed_layer:
            watershed_fid = watershed_feature.GetFID()
            ws_prefix = 'ws_%s_%d' % (watershed_basename, watershed_fid)
            if ws_prefix in scheduled_watershed_prefixes:
                raise ValueError('%s has already been scheduled', ws_prefix)
            scheduled_watershed_prefixes.add(ws_prefix)
            watershed_geom = watershed_feature.GetGeometryRef()
            watershed_area = watershed_geom.GetArea()
            if watershed_area < 0.03:
                #  0.03 square degrees is a healthy underapproximation of
                # 100 sq km which is about the minimum watershed size we'd
                # want.
                continue

            LOGGER.info('processing %s', ws_prefix)
            # make a few subdirectories so we don't explode on number of files per
            # directory. The largest watershed is 726k
            last_digits = '%.4d' % watershed_fid
            local_workspace_dir = os.path.join(
                SDR_WORKSPACES_DIR, last_digits[-1], last_digits[-2],
                last_digits[-3], last_digits[-4],
                "%s" % ws_prefix)
            if not os.path.exists(local_workspace_dir):
                os.makedirs(local_workspace_dir)

            # find EPSG code and pass that/modify SDR for it
            centroid_geom = watershed_geom.Centroid()
            utm_code = (math.floor((centroid_geom.GetX() + 180)/6) % 60) + 1
            lat_code = 6 if centroid_geom.GetY() > 0 else 7
            epsg_code = int('32%d%02d' % (lat_code, utm_code))

            local_watershed_vector_path = os.path.join(
                local_workspace_dir, '%s.gpkg' % ws_prefix)
            make_local_watershed_task = task_graph.add_task(
                func=make_local_watershed,
                args=(
                    watershed_path, watershed_fid, epsg_code,
                    local_watershed_vector_path),
                target_path_list=[local_watershed_vector_path],
                task_name='make local watershed for %s' % ws_prefix)

            # clip dem
            clipped_dir = os.path.join(local_workspace_dir, 'pre_clipped')
            try:
                os.makedirs(clipped_dir)
            except OSError:
                pass
            target_raster_path_list = [
                os.path.join(clipped_dir, '%s_clipped%s.tif' % (
                    raster_type, ws_prefix))
                for raster_type in ['dem', 'erosivity', 'erodibility', 'lulc']]
            base_raster_path_list = [
                dem_vrt_path, erosivity_path, erodibility_path, lulc_path]
            dem_info = pygeoprocessing.get_raster_info(dem_vrt_path)

            dem_pixel_size = dem_info['pixel_size']
            pre_align_task = task_graph.add_task(
                func=pygeoprocessing.align_and_resize_raster_stack,
                args=(
                    base_raster_path_list, target_raster_path_list,
                    ['near']*len(base_raster_path_list),
                    dem_pixel_size, 'intersection'),
                kwargs={
                    'base_vector_path_list': [local_watershed_vector_path],
                    'target_sr_wkt': dem_info['projection']
                    },
                dependent_task_list=[
                    fetch_lulc_task,
                    fetch_erosivity_task,
                    fetch_erodibility_task,
                    make_dem_task,
                    make_local_watershed_task],
                target_path_list=target_raster_path_list,
                task_name='pre-clip for %s' % ws_prefix)

            m_per_deg = length_of_degree(centroid_geom.GetY())
            target_pixel_size = (
                m_per_deg*dem_pixel_size[0], m_per_deg*dem_pixel_size[1])

            sdr_args = {
                'workspace_dir': local_workspace_dir,
                'results_suffix': ws_prefix,
                'dem_path': target_raster_path_list[0],
                'erosivity_path': target_raster_path_list[1],
                'erodibility_path': target_raster_path_list[2],
                'lulc_path': target_raster_path_list[3],
                'watersheds_path': local_watershed_vector_path,
                'biophysical_table_path': biophysical_table_path,
                'threshold_flow_accumulation': 1000,
                'biophysical_table_lucode_header_id': 'ID',
                'k_param': '2',
                'sdr_max': '0.8',
                'ic_0_param': '0.5',
                'local_projection_epsg': epsg_code,
                'target_pixel_size': target_pixel_size,
                'biophysical_table_lucode_field': 'id',
            }
            LOGGER.debug('adding %s', ws_prefix)
            task_graph.add_task(
                func=natcap.invest.sdr.execute,
                args=(sdr_args,),
                target_path_list=[os.path.join(
                    local_workspace_dir,
                    'watershed_results_sdr_%s.shp' % ws_prefix)],
                dependent_task_list=[pre_align_task],
                task_name='sdr for %s' % ws_prefix)

    task_graph.close()
    task_graph.join()


def mosaic_watersheds(task_graph, workspace_dir, base_raster_dir, raster_name):
    """Entry point."""
    try:
        os.makedirs(workspace_dir)
    except OSError:
        pass

    LOGGER.debug("gathering directory list from %s", base_raster_dir)
    leaf_directory_list = (
        (dirpath, filenames) for (dirpath, dirnames, filenames) in os.walk(
            base_raster_dir) if 'intermediate_outputs' in dirnames)

    raster_pattern_to_aggregate = '%s_.*\.tif' % raster_name
    # peek at first element
    sample_dirpath, sample_filenames = next(leaf_directory_list)

    try:
        base_raster_path = next(iter(
            (os.path.join(sample_dirpath, filename)
             for filename in sample_filenames
             if re.match(raster_pattern_to_aggregate, filename))))
    except StopIteration:
        raise ValueError(
            "Expected to find %s in %s but not found" % (
                raster_pattern_to_aggregate, sample_dirpath))

    base_raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    target_raster_path = os.path.join(workspace_dir, '%s.tif' % raster_name)
    target_token_complete_path = '%s_%s.TOKEN' % (
        os.path.splitext(target_raster_path)[0], MOSAIC_DEGREE_CELL_SIZE)
    LOGGER.debug(target_raster_path)
    make_empty_raster_task = task_graph.add_task(
        func=make_empty_wgs84_raster,
        args=(
            MOSAIC_DEGREE_CELL_SIZE, base_raster_info['nodata'][0],
            base_raster_info['datatype'], target_raster_path,
            target_token_complete_path),
        ignore_path_list=[target_raster_path],
        target_path_list=[target_token_complete_path],
        task_name='create empty global %s' % raster_name)
    LOGGER.info("found all the raster suffixes in %s", sample_dirpath)

    leaf_directory_list = (
        (dirpath, filenames) for (dirpath, dirnames, filenames) in os.walk(
            base_raster_dir) if 'intermediate_outputs' in dirnames)
    for dirpath, filenames in leaf_directory_list:
        try:
            base_raster_path = next(iter(
                (os.path.join(dirpath, filename)
                 for filename in filenames
                 if re.match(raster_pattern_to_aggregate, filename))))
        except StopIteration:
            raise RuntimeError(
                "Expected to find %s in %s but not found %s" % (
                    raster_pattern_to_aggregate, dirpath, (
                        dirpath, filenames)))

        target_wgs84_raster_path = '%s_wgs84.tif' % os.path.splitext(
            base_raster_path)[0]
        wgs84_project_task = task_graph.add_task(
            func=pygeoprocessing.warp_raster,
            args=(
                base_raster_path,
                (MOSAIC_DEGREE_CELL_SIZE, -MOSAIC_DEGREE_CELL_SIZE),
                target_wgs84_raster_path, 'near'),
            kwargs={'target_sr_wkt': WSGS84_WKT},
            target_path_list=[target_wgs84_raster_path],
            dependent_task_list=[make_empty_raster_task],
            task_name='wgs84 project %s' % os.path.basename(
                base_raster_path))

        mosaic_complete_token_path = '%s.MOSAICKED' % (
            os.path.splitext(target_wgs84_raster_path)[0])
        mosiac_task = task_graph.add_task(
            func=mosaic_base_into_target,
            args=(
                target_wgs84_raster_path, target_raster_path,
                mosaic_complete_token_path),
            ignore_path_list=[target_raster_path],
            target_path_list=[mosaic_complete_token_path],
            dependent_task_list=[wgs84_project_task],
            task_name='mosiac %s' % (
                os.path.basename(target_wgs84_raster_path)))
        # this ensures that a mosiac will happen one at a time

    target_compressed_path = '%s_compressed.tif' % os.path.splitext(
        target_raster_path)[0]
    task_graph.add_task(
        func=compress_to,
        args=(target_raster_path, 'near', target_compressed_path),
        target_path_list=[target_compressed_path],
        dependent_task_list=[mosiac_task],
        task_name='compress %s' % target_raster_path)

    task_graph.join()
    task_graph.close()


def download_validate_and_unzip(url, target_dir, token_file):
    """Download url to target and write a token file when it unzips."""
    target_path = os.path.join(target_dir, os.path.basename(url))
    url_fetch_and_validate(url, target_path)
    with zipfile.ZipFile(target_path, 'r') as zip_ref:
        zip_ref.extractall(os.path.dirname(target_path))
    with open(token_file, 'w') as token_file:
        token_file.write(str(datetime.datetime.now()))


def url_fetch_and_validate(url, target_path):
    """Download a Google Blob to a given path and hash.

    Parameters:
        url (string): url to a file to fetch.
        target_path (string): path to download the file into, must match
            an embedded hash/algorithm pair.

    Raises:
        ValueError if downloaded file does not match its embedded fingerprint
            where the filename is of the form
            [filename]_[hash_alg]_[fingerprint].ext

    Returns:
        None.

    """
    url_fetcher(url, target_path)
    if not valid_hash(target_path, 'embedded'):
        raise ValueError("%s does not match its expected hash" % target_path)


def url_fetcher(url, path):
    """Download `url` to `path`."""
    LOGGER.info('fetching %s' % path)
    response = urllib.request.urlopen(url)
    with open(path, 'wb') as out_file:
        shutil.copyfileobj(response, out_file)
    response.close()


def valid_hash(file_path, expected_hash, buf_size=2**20):
    """Validate that the file at `file_path` matches `expected_hash`.

    Parameters:
        file_path (str): path to file location on disk.
        expected_hash (str or tuple): if a tuple, a "hash algorithm",
            "expected_hash" pair that will be used to hash `expected_path`
            and confirm that the hash of that file is equivalent to the
            expected hash value. If the `expected_path` does not
            match the hash, this function will raise an AssertionError.

            Otherwise must be the value "embedded" which attempts to parse
            `file_path` for the pattern
            filename_{hash_algorthm}_{hash_value}.{rest of filename}.
        buf_size (int): (optional) number of bytes to read from `file_path`
            at a time for digesting.

    Returns:
        True if `file_path` hashes to `expected_hash`.

    Raises:
        ValueError if `expected_hash == 'embedded'` and `file_path` does not
        match the appropriate file pattern.

        IOError if `file_path` not found.

    """
    if not os.path.exists(file_path):
        raise IOError('%s not found.' % file_path)
    if isinstance(expected_hash, tuple):
        hash_algorithm = expected_hash[0]
        expected_hash_value = expected_hash[1]
    elif expected_hash == 'embedded':
        hash_re_pattern = r'.*_([^_]+)_([0-9a-f]+)\.[^_]*$'
        hash_match = re.match(hash_re_pattern, file_path)
        if not hash_match:
            raise ValueError(
                "file_path: %s did not end "
                "in an [hash_alg]_[hexhash][.ext] format" % {file_path})
        hash_algorithm = hash_match.group(1)
        expected_hash_value = hash_match.group(2)
    else:
        raise ValueError(
            "Invalid value for `expected_hash`, expecting either a tuple "
            "or 'embedded': {expected_hash}")

    try:
        actual_hash = hash_file(file_path, hash_algorithm)
        return expected_hash_value == actual_hash
    except ValueError:
        LOGGER.exception('possible unsupported hash value, not checking')
        return True


def hash_file(file_path, hash_algorithm, buf_size=2**20):
    """Return a hex  digest of `file_path`.

    Parameters:
        file_path (string): path to file to hash.
        hash_algorithm (string): a hash function id that exists in
            hashlib.algorithms_available.
        buf_size (int): number of bytes to read from `file_path` at a time
            for digesting.

    Returns:
        a (hash, crc32) hex digest tuple with hash algorithm `hash_algorithm`
        of the binary contents of `file_path` and the crc32 checksum of that
        file.

    """
    hash_func = hashlib.new(hash_algorithm)
    with open(file_path, 'rb') as f:
        binary_data = f.read(buf_size)
        while binary_data:
            hash_func.update(binary_data)
            binary_data = f.read(buf_size)
    return hash_func.hexdigest()[:32]


def make_vrt(
        base_raster_pattern, target_nodata, target_raster_path,
        target_dem_vrt_token_path):
    """Make a VRT given a list of files.

    Parameters:
        base_raster_pattern (str): pattern of rasters to build to vrt.
        target_nodata (numeric): desired nodata.
        target_raster_path (str): path to desired target vrt
        target_dem_vrt_token_path (str): path to a file to write when
            complete.

    Returns:
        None.

    """
    base_raster_path_list = glob.glob(base_raster_pattern)
    vrt_options = gdal.BuildVRTOptions(VRTNodata=target_nodata)
    gdal.BuildVRT(
        target_raster_path, base_raster_path_list, options=vrt_options)
    target_dem = gdal.OpenEx(target_raster_path, gdal.OF_RASTER)
    if target_dem is not None:
        with open(target_dem_vrt_token_path, 'w') as token_file:
            token_file.write(str(datetime.datetime.now()))
    else:
        raise RuntimeError(
            "didn't make VRT at %s on: %s, from: %s" % (
                target_raster_path, base_raster_path_list,
                base_raster_pattern))


def length_of_degree(lat):
    """Calculate the length of a degree in meters."""
    m1 = 111132.92
    m2 = -559.82
    m3 = 1.175
    m4 = -0.0023
    p1 = 111412.84
    p2 = -93.5
    p3 = 0.118
    lat_rad = lat * math.pi / 180
    latlen = (
        m1 + m2 * math.cos(2 * lat_rad) + m3 * math.cos(4 * lat_rad) +
        m4 * math.cos(6 * lat_rad))
    longlen = abs(
        p1 * math.cos(lat_rad) + p2 * math.cos(3 * lat_rad) + p3 * math.cos(5 * lat_rad))
    return max(latlen, longlen)


def make_empty_wgs84_raster(
        cell_size, nodata_value, target_datatype, target_raster_path,
        target_token_complete_path):
    """Make a big empty raster in WGS84 projection.

    Parameters:
        cell_size (float): this is the desired cell size in WSG84 degree
            units.
        nodata_value (float): desired nodata avlue of target raster
        target_datatype (gdal enumerated type): desired target datatype.
        target_raster_path (str): this is the target raster that will cover
            [-180, 180), [90, -90) with cell size units with y direction being
            negative.
        target_token_complete_path (str): this file is created if the
            mosaic to target is successful. Useful for taskgraph task
            scheduling.

    Returns:
        None.

    """
    gtiff_driver = gdal.GetDriverByName('GTiff')
    try:
        os.makedirs(os.path.dirname(target_raster_path))
    except OSError:
        pass

    n_cols = int(360.0 / cell_size)
    n_rows = int(180.0 / cell_size)

    geotransform = (-180.0, cell_size, 0.0, 90.0, 0, -cell_size)

    target_raster = gtiff_driver.Create(
        target_raster_path, n_cols, n_rows, 1, target_datatype,
        options=(
            'TILED=YES', 'BIGTIFF=YES', 'BLOCKXSIZE=256', 'BLOCKYSIZE=256'))
    target_raster.SetProjection(WSGS84_WKT)
    target_raster.SetGeoTransform(geotransform)
    target_band = target_raster.GetRasterBand(1)
    target_band.SetNoDataValue(nodata_value)
    LOGGER.debug("filling %s with %s" % (target_raster_path, nodata_value))
    target_band.Fill(nodata_value)
    target_band.FlushCache()
    target_band = None
    target_raster = None

    target_raster = gdal.OpenEx(target_raster_path, gdal.OF_RASTER)
    if target_raster:
        with open(target_token_complete_path, 'w') as target_token_file:
            target_token_file.write('complete!')


def make_local_watershed(
        watershed_path, watershed_fid, epsg_code,
        local_watershed_vector_path):
    """Make a local reference for a watershed.

    Parameters:
        watershed_path (str): path to base watershed vector
        watershed_fid (int): FID for the watershed to remove
        epsg_code (str): EPSG code to project local watershed to
        local_watershed_vector_path (str): path to watershed shapefile to
            create.

    Returns:
        None.
    """

    watershed_vector = gdal.OpenEx(watershed_path, gdal.OF_VECTOR)
    watershed_layer = watershed_vector.GetLayer()
    watershed_feature = watershed_layer.GetFeature(watershed_fid)
    watershed_geom = watershed_feature.GetGeometryRef()

    epsg_srs = osr.SpatialReference()
    epsg_srs.ImportFromEPSG(epsg_code)

    wgs84_srs = osr.SpatialReference()
    wgs84_srs.ImportFromEPSG(4326)
    wgs84_to_utm = osr.CoordinateTransformation(wgs84_srs, epsg_srs)

    # clip out watershed to its own file
    # create a new shapefile
    if os.path.exists(local_watershed_vector_path):
        os.remove(local_watershed_vector_path)
    driver = ogr.GetDriverByName('GPKG')
    local_watershed_vector = driver.CreateDataSource(
        local_watershed_vector_path)
    local_watershed_layer = local_watershed_vector.CreateLayer(
        os.path.splitext(os.path.basename(local_watershed_vector_path))[0],
        epsg_srs, ogr.wkbPolygon)
    local_watershed_layer.CreateField(
        ogr.FieldDefn('ws_id', ogr.OFTInteger))
    layer_defn = local_watershed_layer.GetLayerDefn()
    feature_geometry = watershed_geom.Clone()
    watershed_feature = ogr.Feature(layer_defn)
    feature_geometry.Transform(wgs84_to_utm)
    watershed_feature.SetGeometry(feature_geometry)
    watershed_feature.SetField('ws_id', 0)
    local_watershed_layer.CreateFeature(watershed_feature)
    local_watershed_layer.SyncToDisk()
    watershed_geom = None
    feature_geometry = None
    watershed_feature = None
    local_watershed_layer = None
    local_watershed_vector = None


def mosaic_base_into_target(
        base_raster_path, target_raster_path, target_token_complete_path):
    """Copy valid parts of base to target w/r/t correct georeference.

    Parameters:
        base_raster_path (str): a raster with the same cell size,
            coordinate system, and nodata as `target_raster_path`.
        target_raster_path (str): a raster that already exists on disk that
            after this call will contain the non-nodata parts of
            `base_raster_path` that geographically overlap with the target.
        target_token_complete_path (str): this file is created if the
            mosaic to target is successful. Useful for taskgraph task
            scheduling.

    Returns:
        None.

    """
    target_raster = gdal.OpenEx(
        target_raster_path, gdal.OF_RASTER | gdal.GA_Update)
    target_band = target_raster.GetRasterBand(1)
    target_raster_info = pygeoprocessing.get_raster_info(target_raster_path)
    target_nodata = target_raster_info['nodata'][0]
    base_raster_info = pygeoprocessing.get_raster_info(base_raster_path)
    target_gt = target_raster_info['geotransform']
    base_gt = base_raster_info['geotransform']

    target_x_off = int((base_gt[0] - target_gt[0]) / target_gt[1])
    target_y_off = int((base_gt[3] - target_gt[3]) / target_gt[5])

    for offset_dict, band_data in pygeoprocessing.iterblocks(
            (base_raster_path, 1)):
        target_block = target_band.ReadAsArray(
            xoff=offset_dict['xoff']+target_x_off,
            yoff=offset_dict['yoff']+target_y_off,
            win_xsize=offset_dict['win_xsize'],
            win_ysize=offset_dict['win_ysize'])
        valid_mask = numpy.isclose(target_block, target_nodata)
        target_block[valid_mask] = band_data[valid_mask]
        target_band.WriteArray(
            target_block,
            xoff=offset_dict['xoff']+target_x_off,
            yoff=offset_dict['yoff']+target_y_off)
    target_band.FlushCache()
    target_band = None
    target_raster = None

    with open(target_token_complete_path, 'w') as token_file:
        token_file.write('complete!')


def compress_to(base_raster_path, resample_method, target_path):
    """Compress base to target using resample method for overviews."""
    gtiff_driver = gdal.GetDriverByName('GTiff')
    base_raster = gdal.OpenEx(base_raster_path, gdal.OF_RASTER)
    LOGGER.info('compress %s to %s' % (base_raster_path, target_path))
    gtiff_driver.CreateCopy(
        target_path, base_raster, options=(
            'TILED=YES', 'BIGTIFF=YES', 'COMPRESS=LZW',
            'BLOCKXSIZE=256', 'BLOCKYSIZE=256'))
    base_raster = None
    min_dimension = min(
        pygeoprocessing.get_raster_info(target_path)['raster_size'])
    LOGGER.info("min min_dimension %s" % min_dimension)
    raster_copy = gdal.OpenEx(target_path, gdal.OF_RASTER)

    overview_levels = []
    current_level = 2
    while True:
        if min_dimension // current_level == 0:
            break
        overview_levels.append(current_level)
        current_level *= 2
    LOGGER.info('level list: %s' % overview_levels)
    gdal.SetConfigOption('COMPRESS_OVERVIEW', 'LZW')
    raster_copy.BuildOverviews(
        resample_method, overview_levels, callback=_make_logger_callback(
            'build overview for ' + os.path.basename(target_path) +
            '%.2f%% complete'))


if __name__ == '__main__':
    main()

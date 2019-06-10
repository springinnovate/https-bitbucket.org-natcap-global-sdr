"""Script to run a global SDR analysis."""
import multiprocessing
import math
import glob
import urllib
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

DEM_TARGET_NODATA = -32768

N_CPUS = multiprocessing.cpu_count()
TASKGRAPH_REPORTING_FREQUENCY = 5.0
LOGGING_LEVEL = logging.INFO

logging.basicConfig(
    level=LOGGING_LEVEL,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

EROSIVITY_URL = r'https://storage.googleapis.com/global-invest-sdr-data/erosivity_CIAT_50km_md5_8e0d84d5736d118e111b8ee0ded65358.tif'
ERODIBILITY_URL = r'https://storage.googleapis.com/global-invest-sdr-data/erodibility_globe_ISRIC_30arcseconds_md5_e3f8961b77539b686deb9a3d04ee4ce3.tif'
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
        ignore_path_list=[dem_vrt_path],
        target_path_list=[dem_vrt_token_path],
        task_name='make dem vrt')

    task_graph.join()
    scheduled_watershed_prefixes = set()
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
            epsg_srs = osr.SpatialReference()
            epsg_srs.ImportFromEPSG(epsg_code)

            wgs84_srs = osr.SpatialReference()
            wgs84_srs.ImportFromEPSG(4326)
            wgs84_to_utm = osr.CoordinateTransformation(wgs84_srs, epsg_srs)

            # clip out watershed to its own file
            # create a new shapefile
            watershed_vector_path = os.path.join(
                local_workspace_dir, '%s.gpkg' % ws_prefix)
            if os.path.exists(watershed_vector_path):
                os.remove(watershed_vector_path)
            driver = ogr.GetDriverByName('GPKG')
            watershed_vector = driver.CreateDataSource(watershed_vector_path)
            watershed_layer = watershed_vector.CreateLayer(
                os.path.splitext(os.path.basename(watershed_vector_path))[0],
                epsg_srs, ogr.wkbPolygon)
            watershed_layer.CreateField(
                ogr.FieldDefn('ws_id', ogr.OFTInteger))
            layer_defn = watershed_layer.GetLayerDefn()
            feature_geometry = watershed_geom.Clone()
            watershed_feature = ogr.Feature(layer_defn)
            feature_geometry.Transform(wgs84_to_utm)
            watershed_feature.SetGeometry(feature_geometry)
            watershed_feature.SetField('ws_id', 0)
            watershed_layer.CreateFeature(watershed_feature)
            watershed_layer.SyncToDisk()
            watershed_geom = None
            feature_geometry = None
            watershed_feature = None
            watershed_layer = None
            watershed_vector = None

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
                    'base_vector_path_list': [watershed_vector_path],
                    'target_sr_wkt': dem_info['projection']
                    },
                dependent_task_list=[make_dem_task],
                target_path_list=target_raster_path_list,
                task_name='pre-clip for %s' % ws_prefix)

            # clip erosivity
            # clip erodibility
            # cilp lulc

            m_per_deg = length_of_degree(centroid_geom.GetY())
            target_pixel_size = (
                m_per_deg*dem_pixel_size[0], m_per_deg*dem_pixel_size[1])

            # call SDR?
            sdr_args = {
                'workspace_dir': local_workspace_dir,
                'results_suffix': ws_prefix,
                'dem_path': target_raster_path_list[0],
                'erosivity_path': target_raster_path_list[1],
                'erodibility_path': target_raster_path_list[2],
                'lulc_path': target_raster_path_list[3],
                'watersheds_path': watershed_vector_path,
                'biophysical_table_path': biophysical_table_path,
                'threshold_flow_accumulation': 1000,
                'k_param': '2',
                'sdr_max': '0.8',
                'ic_0_param': '0.5',
                'local_projection_epsg': epsg_code,
                'target_pixel_size': target_pixel_size,
                'biophysical_table_lucode_field': 'id',
            }
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
    response = urllib.urlopen(url)
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
            "didn't make VRT at %s on: %s" % (
                target_raster_path, base_raster_path_list))


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


if __name__ == '__main__':
    main()

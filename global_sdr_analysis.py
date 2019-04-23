"""Script to run a global SDR analysis."""
import urllib
import shutil
import re
import hashlib
import datetime
import os
import zipfile
import logging
import sys

from osgeo import gdal
import natcap.invest.sdr
import taskgraph


WORKSPACE_DIR = 'workspace'
CHURN_DIR = os.path.join(WORKSPACE_DIR, 'churn')

N_CPUS = -1
TASKGRAPH_REPORTING_FREQUENCY = 5.0

logging.basicConfig(
    level=logging.DEBUG,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(funcName)s:%(lineno)d] %(message)s'),
    stream=sys.stdout)
LOGGER = logging.getLogger(__name__)

LULC_URL = r''

def main():
    """Entry point."""
    for dir_path in [
            WORKSPACE_DIR, CHURN_DIR]:
        try:
            os.makedirs(dir_path)
        except OSError:
            pass

    task_graph = taskgraph.TaskGraph(
        os.path.join(WORKSPACE_DIR, 'taskgraph_cache'), N_CPUS,
        TASKGRAPH_REPORTING_FREQUENCY)

    lulc_token_path = os.path.join(f'{os.path.basename(LULC_URL)}.COMPLETE')
    fetch_lulc_task = task_graph.add_task(
        n_retries=5,
        func=download_validate_and_unzip,
        args=(LULC_URL, CHURN_DIR, lulc_token_path),
        target_path_list=[lulc_token_path],
        task_name='fetch esacci landuse')

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

    actual_hash = hash_file(file_path, hash_algorithm)
    return expected_hash_value == actual_hash


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


def add_nodata_value(base_raster_path, nodata_value):
    """Set the given nodata value to the raster."""
    raster = gdal.OpenEx(base_raster_path, gdal.OF_RASTER | gdal.GA_Update)
    band = raster.GetRasterBand(1)
    band.SetNoDataValue(nodata_value)
    band.FlushCache()
    band = None
    raster = None


if __name__ == '__main__':
    main()

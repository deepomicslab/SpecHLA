"""
Centralized path resolution for SpecHLA.

Priority chain:
  1. SPECHLA_DB env var (user override)
  2. $CONDA_PREFIX/share/spechla/db (conda install)
  3. Relative to this file (development install)

Usage:
  from spechla_paths import get_db_dir, get_script_dir
"""

import os

_dev_path_done = False

def _setup_dev_path():
    """In dev mode, add repo's bin/ to PATH so vendored binaries are found."""
    global _dev_path_done
    if _dev_path_done:
        return
    _dev_path_done = True
    conda = os.environ.get('CONDA_PREFIX', '')
    if conda and os.path.isdir(os.path.join(conda, 'share', 'spechla')):
        return
    bin_dir = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'bin'))
    if not os.path.isdir(bin_dir):
        return
    path_dirs = os.environ.get('PATH', '').split(os.pathsep)
    if bin_dir in path_dirs:
        return
    dirs = [bin_dir,
            os.path.join(bin_dir, 'SpecHap', 'build'),
            os.path.join(bin_dir, 'extractHairs', 'build'),
            os.path.join(bin_dir, 'fermikit', 'fermi.kit')]
    os.environ['PATH'] = os.pathsep.join(dirs) + os.pathsep + os.environ.get('PATH', '')

def get_db_dir():
    _setup_dev_path()
    env = os.environ.get('SPECHLA_DB', '')
    if env and os.path.isdir(env):
        return env
    conda = os.environ.get('CONDA_PREFIX', '')
    conda_db = os.path.join(conda, 'share', 'spechla', 'db')
    if conda and os.path.isdir(conda_db):
        return conda_db
    return os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'db'))

def get_script_dir():
    _setup_dev_path()
    env = os.environ.get('SPECHLA_SCRIPT', '')
    if env and os.path.isdir(env):
        return env
    conda = os.environ.get('CONDA_PREFIX', '')
    conda_script = os.path.join(conda, 'share', 'spechla', 'script')
    if conda and os.path.isdir(conda_script):
        return conda_script
    return os.path.dirname(os.path.abspath(__file__))

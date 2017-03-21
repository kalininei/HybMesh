# Data will be altered during setup.py procedure and
# written to config_installed.py
# This file is used only for debug purposes
import os.path

exedir = os.path.dirname(__file__)
libdir = os.path.join(exedir, '../../../build/bin')
version = '9.9.9'
last_ver_url = 'https://api.github.com/repos/kalininei/HybMesh/releases/latest'
proj_url = 'https://api.github.com/repos/kalininei/HybMesh/'

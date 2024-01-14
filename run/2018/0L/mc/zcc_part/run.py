#!/usr/bin/env python3

import os

basedir = os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', '..')

import sys
sys.path.insert(0, os.path.join(basedir, 'src'))
from config import Config
from run import run

config = Config('config.yaml')
run(config)

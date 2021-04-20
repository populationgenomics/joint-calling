#!/usr/bin/env python

""" Runs a script as a process
"""

import sys
import subprocess

subprocess.run(sys.argv[1:], check=False)

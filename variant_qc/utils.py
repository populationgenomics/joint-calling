import logging
import os
import subprocess


def run(cmd, silent=False):
    """Run the provided command, logging details and checking for errors.
    """
    if not silent:
        logging.warning(' '.join(str(x) for x in cmd) if not isinstance(cmd, str) else cmd)
    subprocess.check_call(cmd, shell=True, executable=find_bash())


def find_bash():
    for test_bash in ["/bin/bash", "/usr/bin/bash", "/usr/local/bin/bash"]:
        if test_bash and os.path.exists(test_bash):
            return test_bash
    raise IOError("Could not find bash in any standard location. Needed for unix pipes")


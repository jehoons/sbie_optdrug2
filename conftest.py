# -*- coding: utf-8 -*-
#!/usr/bin/python
#*************************************************************************
# Author: {Je-Hoon Song, <song.je-hoon@kaist.ac.kr>
#
# This file is part of {sbie_optdrug}.
#*************************************************************************

import pytest

def pytest_addoption(parser):
    parser.addoption("--with_small", action="store_true",
            help="Use small-scale data for quick test")

    parser.addoption("--force", action="store_true",
            help="Force overwriting existing results")

    parser.addoption("--progress", action="store_true",
            help="display progressbar")


def pytest_generate_tests(metafunc):
    if 'with_small' in metafunc.fixturenames:
        if metafunc.config.option.with_small:
            with_small = True
        else:
            with_small = False

        metafunc.parametrize("with_small", [with_small])

    if 'force' in metafunc.fixturenames:
        if metafunc.config.option.force:
            force = True
        else:
            force = False

        metafunc.parametrize("force", [force])

    if 'progress' in metafunc.fixturenames:
        if metafunc.config.option.force:
            progress = True
        else:
            progress = False

        metafunc.parametrize("progress", [progress])

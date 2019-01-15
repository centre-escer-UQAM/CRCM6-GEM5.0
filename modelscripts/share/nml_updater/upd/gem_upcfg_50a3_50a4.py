#!/usr/bin/env python
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# Author: Stephane Chamberland <stephane.chamberland@canada.ca>
# Copyright: LGPL 2.1
"""
"""

from rpnpy.ftnnml import *

def main(cfg, fromVersion=None, toVersion=None, verbose=False, debug=False):
    """
    """
    if verbose or debug:
        print("Updating from gem 5.0.a3 to 5.0.a4 (%s, %s)" % (fromVersion, toVersion))

    # Special Cases
    # Convert fcst_alarm_s
    # mv: step / fcst_alarm_s > step / step_alarm
    # rm: step / fcst_alarm_s
    dt = None
    try:
        dt = float(cfg.get('step').get('step_dt').get('v').toStr())
    except:
        sys.stdout.write('Error: cannot retrieve dt\n')
    nml1 = 'step'
    try:
        nml = cfg.get(nml1)
    except:
        nml = None
    if nml and dt:
        var1 = 'fcst_alarm_s'
        var2 = 'step_alarm'
        try:
            val = nml.get(var1).get('v').toStr()
            nml.rm(var1)
            units = val.strip().lower()[-1]
            val2 = val.strip().lower()[:-1]
            if units == 's':
                val = float(val2) / float(dt)
            elif units == 'm':
                val = float(val2) * 60. / float(dt)
            elif units == 'h':
                val = float(val2) * 3600. / float(dt)
            elif units == 'd':
                val = float(val2) * 86400. / float(dt)
            else:
                val = float(val)
            nml.add(
                FtnNmlKeyVal(var2, FtnNmlVal(int(val)))
                )
            if verbose or debug:
                sys.stdout.write('Convert: %s/%s => %s\n'
                                 % (nml1, var1, var2))
        except:
            pass

    return cfg

# -*- Mode: C; tab-width: 4; indent-tabs-mode: nil -*-
# vim: set expandtab ts=4 sw=4:
# kate: space-indent on; indent-mode cstyle; indent-width 4; mixedindent off;

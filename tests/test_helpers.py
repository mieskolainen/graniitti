# Helper functions for unit & global tests
#
#
# (c) 2017-2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import subprocess


def execute(cmd, expect="[gr: done]"):
    """
    Execute command
    """
    def exec(c):
        result = subprocess.check_output(c, shell=True, text=True)
        ok = expect in result; assert ok == True
        return result

    if type(cmd) is list:
        for i in range(len(cmd)):
            print(__name__ + f': executing: {cmd[i]}')
            result = exec(cmd[i])
    else:
        result = exec(cmd)

##############################################################################
# Medical Image Registration ToolKit (MIRTK)
#
# Copyright 2017 Imperial College London
# Copyright 2017 Andreas Schuh
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
##############################################################################

"""Auxiliary functions for batch execution of MIRTK commands using SLURM."""

import re
import os
import sys
import subprocess
import time


# ----------------------------------------------------------------------------
def submit(name, command=None, args=[], opts={}, script=None, tasks=0, deps=[],
           logdir=None, log=None, queue='long', threads=8, memory=8 * 1024,
           workdir=None, verbose=1):
    """Submit batch job to SLURM."""
    if threads <= 0:
        raise ValueError("Must specify number of threads when executing as SLURM batch job!")
    if script:
        if command:
            raise ValueError("Keyword arguments 'command' and 'script' are mutually exclusive")
        shexec = "#!/bin/bash\nexec {0} <(cat <<END_OF_SCRIPT\n".format(sys.executable)
        shexec += script.format(**opts)
        shexec += "\nEND_OF_SCRIPT)"
        if tasks > 0:
            shexec += " $SLURM_ARRAY_TASK_ID"
        shexec += "\n"
        script = shexec
    else:
        script = "#!/bin/bash\n"
        script += "\"{0}\" -c 'import sys; import socket;".format(sys.executable)
        script += " sys.stdout.write(\"Host: \" + socket.gethostname() + \"\\n\\n\");"
        script += " sys.path.insert(0, \"{0}\");".format(os.path.dirname(os.path.dirname(__file__)))
        script += " import mirtk; mirtk.check_call([\"{0}\"] + sys.argv[1:])".format(command if command else name)
        script += "'"
        for arg in args:
            arg = str(arg)
            if ' ' in arg:
                arg = '"' + arg + '"'
            script += ' ' + arg
        for opt in opts:
            arg = opts[opt]
            if opt[0] != '-':
                opt = '-' + opt
            script += ' ' + opt
            if arg is not None:
                if isinstance(arg, (list, tuple)):
                    arg = ' '.join([str(x) for x in arg])
                else:
                    arg = str(arg)
                script += ' ' + arg
        script += " -threads {0}".format(threads)
        script += "\n"
    argv = [
        'sbatch',
        '-J', name,
        '-n', '1',
        '-c', str(threads),
        '-p', queue,
        '--mem={0}M'.format(memory)
    ]
    if tasks > 0:
        argv.append('--array=0-{}'.format(tasks - 1))
    if logdir or log:
        if not logdir:
            logdir = os.path.dirname(log)
        elif not log:
            log = os.path.join(logdir, name)
            if tasks > 0:
                log += "_%A.%a.log"
            log += "_%j.log"
        if not os.path.exists(logdir):
            os.makedirs(logdir)
        argv.extend(['-o', log, '-e', log])
    if deps:
        if isinstance(deps, int):
            deps = [deps]
        deps = [str(dep) for dep in deps if dep > 0]
        if deps:
            argv.append('--dependency=afterok:' + ',afterok:'.join(deps))
    if workdir:
        argv.append('--workdir=' + os.path.abspath(workdir))
    proc = subprocess.Popen(argv, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    (out, err) = proc.communicate(input=script.encode('utf-8'))
    if proc.returncode != 0:
        raise Exception(err)
    match = re.match('Submitted batch job ([0-9]+)', out)
    if not match:
        raise Exception("Failed to determine job ID from sbatch output:\n" + out)
    jobid = int(match.group(1))
    if verbose > 0:
        if tasks > 0:
            print("  Submitted batch {} (JobId={}, Tasks={})".format(name, jobid, tasks))
        else:
            print("  Submitted job {} (JobId={})".format(name, jobid))
    return jobid


# ----------------------------------------------------------------------------
def wait(jobs, max_time=0, interval=60, verbose=0):
    """Wait for SLURM jobs to finish."""
    if not isinstance(jobs, (list, tuple)):
        jobs = [jobs]
    jobs = [job for job in jobs if job > 0]
    num_wait = len(jobs)
    if num_wait == 0:
        return True
    num_fail = 0
    total_time = 0
    re_state = re.compile("JobState=([A-Z]+)")
    jobs = [str(job) for job in jobs]
    iterations = 0
    while num_wait > 0 and (max_time <= 0 or total_time < max_time):
        time.sleep(interval)
        total_time += interval
        done = []
        fail = []
        batch = []
        table = subprocess.check_output(["sacct", "--brief", "--parsable", "--jobs={}".format(','.join(jobs))])
        for line in table.splitlines():
            cols = line.split("|")
            if cols[0] in jobs:
                if cols[1] == "COMPLETED":
                    done.append(cols[0])
                elif cols[1] not in ("PENDING", "SUSPENDED", "RUNNING"):
                    fail.append(cols[0])
            elif cols[0].endswith(".batch"):
                batch.append(cols[0][0:-6])
        num_jobs = len(jobs)
        num_done = 0
        num_fail = 0
        for job in batch:
            try:
                info = subprocess.check_output(["scontrol", "show", "job", job])
                for line in info.splitlines():
                    match = re_state.search(line)
                    if match:
                        num_jobs += 1
                        if match.group(1) == "COMPLETED":
                            num_done += 1
                        elif match.group(1) not in ("PENDING", "SUSPENDED", "RUNNING"):
                            num_fail += 1
                num_jobs -= 1
                try:
                    done.remove(job)
                except ValueError:
                    pass
                try:
                    fail.remove(job)
                except ValueError:
                    pass
            except subprocess.CalledProcessError:
                pass  # scontrol forgets about no longer queued/running jobs
        num_done += len(done)
        num_fail += len(fail)
        num_wait = num_jobs - num_done - num_fail
        if verbose > 0 and (num_wait <= 0 or iterations % verbose == 0):
            sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
            sys.stdout.write(" WAIT {} job(s) running/suspended/pending".format(num_wait))
            if num_fail > 0:
                sys.stdout.write(", {} failed".format(num_fail))
            sys.stdout.write("\n")
            sys.stdout.flush()
        ++iterations
    if num_wait > 0 and max_time > 0 and total_time >= max_time:
        raise Exception("Exceeded maximum time waiting for jobs to complete!")
    if total_time > 0:
        time.sleep(10)  # wait a bit for files to be available from all NFS clients
    return num_fail == 0

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

"""Auxiliary functions for batch execution of MIRTK commands using HTCondor."""

import re
import os
import sys
import subprocess
from xml.etree import ElementTree


# ----------------------------------------------------------------------------
def submit(name, command=None, args=[], opts={}, script=None, tasks=0, deps=[],
           threads=0, memory=8 * 1024, retries=5, requirements=[], environment={},
           logdir=None, log=None, workdir=None, verbose=1):
    """Submit batch job to HTCondor."""
    if deps:
        raise NotImplementedError("Cannot submit individual HTCondor jobs with dependencies yet, this requires use of DAGMan")
    if logdir or log:
        if not logdir:
            log = os.path.abspath(log)
            logdir = os.path.dirname(log)
        elif not log:
            logdir = os.path.abspath(logdir)
            if tasks > 0:
                log = os.path.join(logdir, name + "_$(Cluster).$(Process).log")
            else:
                log = os.path.join(logdir, name + "_$(Cluster).log")
        os.makedirs(logdir)
    jobdesc = "universe = vanilla\n"
    if threads > 0:
        jobdesc += "request_cpus = {0}\n".format(threads)
    if memory > 0:
        jobdesc += "request_memory = {0}\n".format(memory)
    if requirements:
        jobdesc += "requirements = " + " && ".join(requirements) + "\n"
    if environment:
        jobdesc += "environment = \""
        for envname, envval in environment.items():
            jobdesc += " {0}='{1}'".format(envname, ':'.join(envval) if isinstance(envval, (list, tuple)) else envval)
        jobdesc += "\"\n"
    if workdir:
        jobdesc += "initialdir = {0}\n".format(os.path.abspath(workdir))
    # Note: MIRTK executables return exit code 6 when memory allocation fails, other codes are kill/term signals
    jobdesc += "on_exit_remove = (ExitBySignal == False && ExitCode != 6 && ExitCode != 247 && ExitCode != 241) || (ExitBySignal == True && ExitSignal != 9 && ExitSignal != 15)\n"
    if retries > 0:
        jobdesc += "max_retries = {0}\n".format(retries)
    if script:
        if command:
            raise ValueError("Keyword arguments 'command' and 'script' are mutually exclusive")
        if not log:
            raise ValueError("Script submission of batch to HTCondor requires log path for script file!")
        script_path = os.path.join(logdir, name + ".py")
        with open(script_path, "wt") as f:
            f.write(script.format(**opts))
        jobdesc += "executable = {0}\n".format(sys.executable)
        jobdesc += "arguments = \"'{0}'".format(script_path)
        if tasks > 0:
            jobdesc += " $(Process)"
        jobdesc += "\"\n"
        if log:
            jobdesc += "output = {0}\n".format(log)
            jobdesc += "error = {0}\n".format(log)
        jobdesc += "queue"
        if tasks > 0:
            jobdesc += " {}".format(tasks)
        jobdesc += "\n"
    else:
        jobdesc += "executable = {0}\n".format(sys.executable)
        jobdesc += "arguments = \"-c 'import sys; import socket;"
        jobdesc += " sys.stdout.write(\"\"Host: \"\" + socket.gethostname() + \"\"\\n\\n\"\");"
        jobdesc += " sys.path.insert(0, \"\"{0}\"\");".format(os.path.dirname(os.path.dirname(__file__)))
        jobdesc += " import mirtk; mirtk.check_call([\"\"{0}\"\"] + sys.argv[1:])'".format(command if command else name)
        for arg in args:
            arg = str(arg)
            if ' ' in arg:
                arg = "'" + arg + "'"
            jobdesc += ' ' + arg
        for opt in opts:
            arg = opts[opt]
            if opt[0] != '-':
                opt = '-' + opt
            jobdesc += ' ' + opt
            if arg is not None:
                if isinstance(arg, (list, tuple)):
                    arg = ' '.join([str(x) for x in arg])
                else:
                    arg = str(arg)
                jobdesc += ' ' + arg
        jobdesc += "\"\n"
        if log:
            jobdesc += "output = {0}\n".format(log)
            jobdesc += "error = {0}\n".format(log)
        jobdesc += "queue\n"
    proc = subprocess.Popen(["condor_submit", "-"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    (out, err) = proc.communicate(input=jobdesc.encode('utf-8'))
    if proc.returncode != 0:
        raise Exception(err)
    match = re.search('[0-9]+ job(s) submitted to cluster ([0-9]+)\.', out)
    if not match:
        match = re.search('\*\* Proc ([0-9]+)(\.[0-9]+)?:', out)
        if not match:
            raise Exception("Failed to determine job ID from condor_submit output:\n" + out)
    jobid = int(match.group(1))
    if verbose > 0:
        if tasks > 0:
            print("  Submitted batch {} (JobId={}, Tasks={})".format(name, jobid, tasks))
        else:
            print("  Submitted job {} (JobId={})".format(name, jobid))
    return jobid if tasks == 0 else (jobid, tasks)


# ----------------------------------------------------------------------------
def _parse_condor_xml(s):
    """Parse XML ClassAds returned by condor_q or condor_history with -xml option."""
    # Note: condor_history may return multiple (only one non-empty) <classads> tags
    end = 0
    classads = None
    starttag = "<classads>"
    endtag = "</classads>"
    while True:
        start = s.find(starttag, end)
        if start == -1:
            break
        end = s.find(endtag, start)
        if end == -1:
            raise ValueError("Malformed <classads> XML, could not find matching </classads>!")
        end += len(endtag)
        elem = ElementTree.fromstring(s[start:end])
        if classads is None:
            classads = elem
        else:
            classads.extend(elem)
    return classads


# ----------------------------------------------------------------------------
def wait(jobs, max_time=0, max_error=5, interval=60, verbose=0):
    if not isinstance(jobs, list):
        jobs = [jobs]
    clusters = []
    for job in jobs:
        if isinstance(job, tuple):
            cluid = job[0]
            tasks = job[1]
        else:
            cluid = job
            tasks = 1
        if cluid > 0:
            clusters.append((cluid, tasks))
    num_wait = len(clusters)
    num_error = 0
    total_time = 0
    iterations = 0
    while num_wait > 0 and (max_time <= 0 or total_time < max_time or total_time == 0) and num_error < max_error:
        time.sleep(interval)
        total_time += interval
        try:
            num_pending = 0
            num_running = 0
            num_suspended = 0
            num_held = 0
            num_done = 0
            for cluster, tasks in clusters:
                classads = subprocess.check_output(["condor_q", "-xml", str(cluster)], stderr=subprocess.STDOUT)
                classads = _parse_condor_xml(classads)
                for process in range(tasks):
                    classad = classads.find(".c/a[@n='ClusterId'][i='{0}']/../a[@n='ProcId'][i='{1}']/..".format(cluster, process))
                    if classad is None:
                        num_done += 1
                    else:
                        status = int(classad.find("a[@n='JobStatus']/i").text)
                        # 1) Idle
                        # 2) Running
                        # 3) Removed
                        # 4) Completed (also when failed, check ExitCode afterwards using condor_history)
                        # 5) Held
                        # 6) Transferring Output
                        # 7) Suspended
                        if status == 1:
                            num_pending += 1
                        elif status == 2 or status == 6:
                            num_running += 1
                        elif status == 3:
                            num_done += 1
                        elif status == 4:
                            num_done += 1
                        elif status == 5:
                            num_held += 1
                        elif status == 7:
                            num_suspended += 1
                        else:
                            raise Exception("Unknown job status: {}".format(status))
            num_wait = num_running + num_pending + num_held + num_suspended
            if verbose > 0 and (num_wait <= 0 or iterations % verbose == 0):
                sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
                sys.stdout.write(" WAIT {p} pending, {r} running, {s} suspended, {h} held, {d} completed\n".format(
                    p=num_pending, r=num_running, s=num_suspended, h=num_held, d=num_done
                ))
                sys.stdout.flush()
            num_error = 0
        except subprocess.CalledProcessError:
            sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
            sys.stdout.write(" WAIT Failed to retrieve job status, will retry {0} more times!\n".format(max_error - num_error))
            sys.stdout.flush()
            num_error += 1
        iterations += 1
    if num_error >= max_error:
        raise Exception("Exceeded maximum number of retries to query status of jobs!")
    if num_wait > 0 and max_time > 0 and total_time >= max_time:
        raise Exception("Exceeded maximum time waiting for jobs to complete!")
    num_fail = 0
    if total_time > 0:
        time.sleep(10)
        num_jobs = 0
        unknown = {}
        for cluster, tasks in clusters:
            num_jobs += tasks
            unknown[cluster] = [0] * tasks
        num_error = 0
        num_unknown = num_jobs
        while num_unknown > 0 and num_error <= max_error:
            try:
                num_fail = 0
                for cluster, tasks in clusters:
                    classads = _parse_condor_xml(subprocess.check_output(["condor_history", "-xml", str(cluster)]))
                    for process in range(tasks):
                        classad = classads.find(".c/a[@n='ClusterId'][i='{0}']/../a[@n='ProcId'][i='{1}']/..".format(cluster, process))
                        if classad is None:
                            unknown[cluster][process] += 1
                        else:
                            unknown[cluster][process] = 0
                            status = int(classad.find("a[@n='JobStatus']/i").text)
                            exit_code = int(classad.find("a[@n='ExitCode']/i").text)
                            if status != 4 or exit_code != 0:
                                num_fail += 1
                num_unknown = 0
                for cluster, tasks in clusters:
                    for process in range(tasks):
                        if unknown[cluster][process] > 0:
                            if unknown[cluster][process] > max_error:
                                raise Exception("Could not retrieve exit code of job {}.{} for the past {} attempts!".format(cluster, process, unknown[cluster][process]))
                            num_unknown += 1
                if verbose > 0:
                    sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
                    sys.stdout.write(" DONE {0} succeeded, {1} failed".format(num_jobs - num_fail - num_unknown, num_fail))
                    if num_unknown > 0:
                        sys.stdout.write(", {} unknown, will retry".format(num_unknown))
                    sys.stdout.write("\n")
                    sys.stdout.flush()
            except subprocess.CalledProcessError:
                sys.stdout.write("{:%Y-%b-%d %H:%M:%S}".format(datetime.now()))
                sys.stdout.write(" WAIT Failed to retrieve exit codes, will retry {0} more times!\n".format(max_error - num_error))
                sys.stdout.flush()
                num_error += 1
        if num_error >= max_error:
            raise Exception("Exceeded maximum number of retries to query exit codes of jobs!")
    return num_fail == 0

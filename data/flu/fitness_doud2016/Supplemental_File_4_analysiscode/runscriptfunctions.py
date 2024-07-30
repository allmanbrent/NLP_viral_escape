import os
import sys
import string
import time
import multiprocessing
import glob

def RunScript(rundir, run_name, script_name, commands, use_sbatch, sbatch_cpus, walltime=None):
    """Runs a script with specified command-line options.

    *rundir* is the directory in which we run the job. Created if it does
    not exist.
    *run_name* is the name of the run, which should be a string without
    spaces.
    *script_name* is the name of the script that we run.
    *commands* is a list of commands for the script (command-line arguments).
    *use_sbatch* is a Boolean switch specifying whether we use ``sbatch``
    to run the script. If *False*, the script is just run with the command
    line instruction. If *True*, then ``sbatch`` is used, and the command file
    has the prefix *run_name* followed by the suffix ``.sbatch``.
    *sbatch_cpus* is an option that is only meaningful if *use_sbatch* is 
    *True*. It gives the integer number of CPUs that are claimed via
    ``sbatch`` using the option ``sbatch -c``. 
    *waltime* is an option that is only meaningful if *use_sbatch* is
    *True*. If so, it should be an integer giving the number of hours 
    to allocate for the job. If *walltime* has its default value of 
    *None*, no wall time for the job is specified.
    It is assumed that the script can be run at the command line using::
        script_name (command-line options)
    Returns *runfailed*: *True* if run failed, and *False* otherwise.
    """


    print "Running %s for %s in directory %s..." % (script_name, run_name, rundir)
    currdir = os.getcwd()
    if not os.path.isdir(rundir):
        os.mkdir(rundir)
    os.chdir(rundir)

    if (not run_name) or not all([x not in string.whitespace for x in run_name]):
        raise ValueError("Invalid run_name of %s" % run_name)

    if use_sbatch:
        sbatchfile = '%s.sbatch' % run_name # sbatch command file
        jobidfile = 'sbatch_%s_jobid' % run_name # holds sbatch job id
        jobstatusfile = 'sbatch_%s_jobstatus' % run_name # holds sbatch job status
        joberrorsfile = 'sbatch_%s_errors' % run_name # holds sbatch job errors
        sbatch_f = open(sbatchfile, 'w')
        sbatch_f.write('#!/bin/sh\n#SBATCH\n')
        if walltime:
            sbatch_f.write('#PBS -l walltime=%d:00:00\n' % walltime)


        sbatch_f.write('%s ' % (script_name))
        for command in commands:
            sbatch_f.write('%s ' % command)
        sbatch_f.close()


        os.system('sbatch -c %d -e %s %s > %s' % (sbatch_cpus, joberrorsfile, sbatchfile, jobidfile))

        time.sleep(10) # short delay
        jobid = int(open(jobidfile).read().split()[-1])
        nslurmfails = 0
        while True:
            time.sleep(30) # delay
            returncode = os.system('squeue -j %d > %s' % (jobid, jobstatusfile))
            if returncode != 0:
                nslurmfails += 1
                if nslurmfails > 180: # error is squeue fails at least 180 consecutive times
                    raise ValueError("squeue is continually failing, which means that slurm is not working on your system. Note that although this script has crashed, many of the jobs submitted via slurm may still be running. You'll want to monitor (squeue) or kill them (scancel) -- unfortunately you can't do that until slurm starts working again.")
                continue # we got an error while trying to run squeue
            nslurmfails = 0
            lines = open(jobstatusfile).readlines()
            if len(lines) < 2:
                break # no longer in slurm queue
        errors = open(joberrorsfile).read().strip()
    
    else: 
        full_command = '%s ' % (script_name) + ' '.join(commands)
        errors = os.system(full_command)

    os.chdir(currdir)
    if errors:
        print "Done running %s for %s in directory %s.\n" % (script_name, run_name, rundir)
    else:
        print "Done running %s for %s in directory %s.\n" % (script_name, run_name, rundir)


def RunProcesses(processes, nmultiruns):
    """Runs a list *multiprocessing.Process* processes.
    *processes* is a list of *multiprocessing.Process* objects that
    have not yet been started.
    *nmultiruns* is an integer >= 1 indicating the number of simultaneous
    processes to run.
    Runs the processes in *processes*, making sure to never have more than
    *nmultiruns* running at a time. If any of the processes fail (return
    an exitcode with a boolean value other than *False*), an exception
    is raised immediately. Otherwise, this function finishes when all
    processes have completed.
    """
    if not (nmultiruns >= 1 and isinstance(nmultiruns, int)):
        raise ValueError("nmultiruns must be an integer >= 1")
    processes_started = [False] * len(processes)
    processes_running = [False] * len(processes)
    processes_finished = [False] * len(processes)
    while not all(processes_finished):
        if (processes_running.count(True) < nmultiruns) and not all(processes_started):
            i = processes_started.index(False)
            processes[i].start()
            processes_started[i] = True
            processes_running[i] = True
        for i in range(len(processes)):
            if processes_running[i]:
                if not processes[i].is_alive():
                    processes_running[i] = False
                    processes_finished[i] = True
                    if processes[i].exitcode:
                        raise IOError("One of the processes failed to complete.")
        time.sleep(45)

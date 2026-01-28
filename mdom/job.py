#!/usr/bin/env python

# Sample script for use on npx-submitter
# Jessie Thwaites, Nov 2025

# Use htcondor2
import htcondor2
# import numpy as np
import argparse, os, pwd

parser =  argparse.ArgumentParser('Basic job example')
parser.add_argument('--wait', default=10, type=float, 
                    help='sleep time to pass (seconds)')
parser.add_argument('--njobs', default=2, type=int,
                    help='number of jobs to run')
parser.add_argument('--no_dag', default=False, action='store_true',
                    help='Submit jobs without a DAG (not recommended for large numbers of jobs)')
parser.add_argument('--cwd', type=str, default=os.getcwd(),
                    help='Directory to work out of. Expects executable to be in this directory')
# currently, the main script that is run is in the parent directory of cwd. you can change this here
parser.add_argument('--script_path', type=str, default=os.path.join(os.path.dirname(os.getcwd()), 'sleep.py'),
                    help='Path to main script run in the job')
# the virtual environment is also expected to be in this parent directory of cwd. can change that here as well
parser.add_argument('--venv_dir', type=str, default=os.path.join(os.path.dirname(os.getcwd()), 'test_venv'),
                    help='Path to virtual environment you want to use (expects main directory and sources bin/activate within this)')
args = parser.parse_args()

# Config of the condor jobs
# remember to make the executable have execute permissions! "chmod a+x <executable>" at the command line if needed before running
username = pwd.getpwuid(os.getuid())[0]
executable = os.path.join(args.cwd, "job.sh")
condor_config = {
        "executable": executable,
        "arguments": "-s $(script) -w $(wait) -j $(njob) -v $(venv)",          # we will pass in the value for this macro via itemdata
        "should_transfer_files": "yes",                 # copy error and output files
        "scratchdir":f"/scratch/{username}/test_job/",  # directory to save logs/out/error files - should always be in your /scratch/user!
        "output": "$(scratchdir)/out/job.$(ClusterId).$(ProcId).out",
        "error": "$(scratchdir)/err/job.$(ClusterId).$(ProcId).err",
        "log": "$(scratchdir)/log/job.$(ClusterId).log",
        "request_cpus": "4",                            # CPU request for each job
        "request_memory":"4GB"                          # memory request for each job. 1GB is the minimum
        # 'accounting_group':'1_week'                   # default is 48 hours max runtime, add to 1 week queue if you need a long job
}

# make all the output/log/err directories
if not os.path.exists(condor_config["scratchdir"]):
    os.makedirs(condor_config["scratchdir"])
if not os.path.exists(os.path.join(condor_config["scratchdir"], 'err')):
    os.mkdir(os.path.join(condor_config["scratchdir"], 'err'))
if not os.path.exists(os.path.join(condor_config["scratchdir"], 'out')):
    os.mkdir(os.path.join(condor_config["scratchdir"], 'out'))
if not os.path.exists(os.path.join(condor_config["scratchdir"], 'log')):
    os.mkdir(os.path.join(condor_config["scratchdir"], 'log'))

print('Job logs, error, and output will be written to: ', condor_config["scratchdir"])
submit_config = htcondor2.Submit(condor_config)

job_args = []
# Create the dictionary of args to be passed
for n in range(args.njobs):
    # wait and cwd are arguments that are the same across all jobs, while we iterate over njob
    # htcondor expects all args to be strings
    # pass the script and the virtual env to the executable
    job_args.append({'script':args.script_path, 'wait':str(args.wait), 'njob': str(n), 'venv':args.venv_dir})

if args.no_dag:

    schedd = htcondor2.Schedd()
    print('Submitting {} jobs.'.format(len(job_args)))
    submit_result = schedd.submit(submit_config, itemdata = iter(job_args))  # submit one job for each arg option
    print('Submitted to cluster {}.'.format(submit_result.cluster()))

else:
    """
    import htcondor2.dags # this has to be imported directly for the dags to work
    import shutil

    dag = htcondor2.dags.DAG()
    dag_dir = os.path.join(condor_config['scratchdir'], 'dag')

    # create the job object inside the node
    # you can also define child layers similarly
    test_dag = dag.layer(
        name = os.path.join(dag_dir,'test_job_dag'), # had to use absolute paths to get this to work
        submit_description = submit_config,
        vars = job_args,
    )

    # all files written by condor must be in /scratch - including your dag
    dag_file = htcondor2.dags.write_dag(dag, dag_dir)
    os.chdir(condor_config['scratchdir'])

    dag_submit = htcondor2.Submit.from_dag(str(dag_file), 
                                           {'force': 1}) # this ignores any rescue dags. comment out to run rescue dag
    # htcondor expects your executable to be right next to your dag
    shutil.copy2(executable, dag_dir) 

    schedd = htcondor2.Schedd()
    print('Submitting {} jobs.'.format(len(job_args)))
    submit_result = schedd.submit(dag_submit)
    print("DAGMan job cluster is {}".format(submit_result.cluster()))
    """
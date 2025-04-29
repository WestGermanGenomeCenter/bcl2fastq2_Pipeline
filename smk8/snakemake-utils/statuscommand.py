#!/usr/bin/env python3
# This is a small script to allow Snakemake to query the PBSPro scheduler
# for job status. Use with --cluster-status "python statuscommand.py".
# Author: Lukas Rose <lukas.rose@uni-duesseldorf.de>

import sys, logging, argparse, subprocess, json
from enum import Enum
from logging import StreamHandler, FileHandler

class JobState(Enum):
	running = 0
	success = 1
	failed = -1

# Source: PBS Professional 18.2, Administrator's Guide, Section 14.9: Job Exit Status Code
PBS_EXIT_CODES = {
	0: {
		"name": "JOB_EXEC_OK",
		"description": "Job execution was successful",
		"state": JobState.success
	},
	-1: {
		"name": "JOB_EXEC_FAIL1",
		"description": "Job execution failed, before files, no retry",
		"state": JobState.failed
	},
	-2: {
		"name": "JOB_EXEC_FAIL2",
		"description": "Job execution failed, after files, no retry",
		"state": JobState.failed
	},
	-3: {
		"name": "JOB_EXEC_RETRY",
		"description": "Job execution failed, do retry",
		"state": JobState.failed
	},
	-4: {
		"name": "JOB_EXEC_INITABT",
		"description": "Job aborted on MoM initialization",
		"state": JobState.failed
	},
	-5: {
		"name": "JOB_EXEC_INITRST",
		"description": "Job aborted on MoM initialization, checkpoint, no migrate",
		"state": JobState.failed
	},
	-6: {
		"name": "JOB_EXEC_INITRMG",
		"description": "Job aborted on MoM initialization, checkpoint, ok migrate",
		"state": JobState.failed
	},
	-7: {
		"name": "JOB_EXEC_BADRESRT",
		"description": "Job restart failed",
		"state": JobState.failed
	},
	-8: {
		"name": "JOB_EXEC_GLOBUS_INIT_RETRY",
		"description": "Globus can still send jobs to PBS, but PBS no longer supports sending jobs to Globus. No longer used. Initialization of Globus job failed; do retry",
		"state": JobState.failed
	},
	-9: {
		"name": "JOB_EXEC_GLOBUS_INIT_FAIL",
		"description": "Globus can still send jobs to PBS, but PBS no longer supports sending jobs to Globus. No longer used. Initialization of Globus job failed; no retry",
		"state": JobState.failed
	},
	-10: {
		"name": "JOB_EXEC_FAILUID",
		"description": "Invalid UID/GID for job",
		"state": JobState.failed
	},
	-11: {
		"name": "JOB_EXEC_RERUN",
		"description": "Job was rerun",
		"state": JobState.failed
	},
	-12: {
		"name": "JOB_EXEC_CHKP",
		"description": "Job was checkpointed and killed",
		"state": JobState.failed
	},
	-13: {
		"name": "JOB_EXEC_FAIL_PASSWORD",
		"description": "Job failed due to a bad password",
		"state": JobState.failed
	},
	-14: {
		"name": "JOB_EXEC_RERUN_ON_SIS_FAIL",
		"description": "Job was requeued (if rerunnable) or deleted (if not) due to a communication failure between Mother Superior and a Sister",
		"state": JobState.failed
	},
	-15: {
		"name": "JOB_EXEC_QUERST",
		"description": "Requeue job for restart from checkpoint",
		"state": JobState.failed
	},
	-16: {
		"name": "JOB_EXEC_FAILHOOK_RERUN",
		"description": "Job execution failed due to hook rejection; requeue for later retry",
		"state": JobState.failed
	},
	-17: {
		"name": "JOB_EXEC_FAILHOOK_DELETE",
		"description": "Job execution failed due to hook rejection; delete the job at end",
		"state": JobState.failed
	},
	-18: {
		"name": "JOB_EXEC_HOOK_RERUN",
		"description": "A hook requested for job to be requeued",
		"state": JobState.failed
	},
	-19: {
		"name": "JOB_EXEC_HOOK_DELETE",
		"description": "A hook requested for job to be deleted",
		"state": JobState.failed
	},
	-20: {
		"name": "JOB_EXEC_RERUN_MS_FAIL",
		"description": "Mother superior connection failed",
		"state": JobState.failed
	},
	"failed": {
		"name": "JOB_SCRIPT_FAILED",
		"description": "The exit value of the jobscript was {exit_code}. Probably one of the commands failed",
		"state": JobState.failed
		
	},
	"killed": {
		"name": "JOB_SCRIPT_KILLED",
		"description": "The job was killed with signal {signal}. See 'kill -l' for a list of signal names on your system",
		"state": JobState.failed
	},
	"default": {
		"name": "EXIT_UNKNOWN",
		"description": "An unknown exit code was encountered: {exit_code}. Ask the system administrator for help",
		"state": JobState.failed
	}
}

# Source: PBS Professional 18.2, Reference Guide, 2.56.3: Options to qselect, Table 2-23: Job States
PBS_JOB_STATES = {
	"B": {
		"name": "STARTED",
		"description": "Job array has started execution",
		"state": JobState.running
	},
	"E": {
		"name": "EXITING",
		"description": "Job is exiting",
		"state": JobState.running
	},
	"H": {
		"name": "HELD",
		"description": "Job is held",
		"state": JobState.running
	},
	"M": {
		"name": "MOVED",
		"description": "Job is moved",
		"state": JobState.running
	},
	"Q": {
		"name": "QUEUED",
		"description": "Job is queued and waiting to start",
		"state": JobState.running
	},
	"R": {
		"name": "RUNNING",
		"description": "Job is currently running",
		"state": JobState.running
	},
	"S": {
		"name": "SUSPENDED",
		"description": "Job is suspended",
		"state": JobState.running
	},
	"T": {
		"name": "TRANSITING",
		"description": "Job is transiting",
		"state": JobState.running
	},
	"U": {
		"name": "WAIT_USER",
		"description": "Job suspended due to workstation user activity",
		"state": JobState.running
	},
	"W": {
		"name": "WAIT",
		"description": "Job is waiting",
		"state": JobState.running
	},
	"X": {
		"name": "EXITED",
		"description": "The eXited state. Subjobs only",
		"state": JobState.running
	},
	"default": {
		"name": "STATUS_UNKNOWN",
		"description": "An unknown job status was encountered: {job_status}. Ask the system administrator for help",
		"state": JobState.running
	},
	"unknown": {
		"name": "JOB_UNKNOWN",
		"description": "An unknown job ID was encountered: {job_id}. Make sure a job with the given ID exists.",
		"state": JobState.failed
	},
	"error": {
		"name": "UNKNOWN_ERROR",
		"description": "An unknown error ocurred while getting the job state for job {job_id}: {ex}. Please contact the system administrator.",
		"state": JobState.failed
	}
}

def query_qstat(job_id):
	job_data = None
	try:
		qstat_result = subprocess.run(["qstat", "-f", "-x", "-F", "json", str(job_id)], stdout=subprocess.PIPE)
		job_data = json.loads(qstat_result.stdout)
	except Exception as ex:
		logging.getLogger().error("An exception occurred when querying qstat: {ex}".format(ex=ex))
	return job_data

def decode_state_dict(state_dict, job_state, **kwargs):
	state = state_dict[job_state]["state"]
	name = state_dict[job_state]["name"]
	description = state_dict[job_state]["description"].format(**kwargs)
	
	return (state, name, description)

def decode_job_status(job_id, job_data):
	try:
		if ("Jobs" in job_data and job_id in job_data["Jobs"]):
			job_state = job_data["Jobs"][job_id]["job_state"]
			
			if(job_state in PBS_JOB_STATES):
				result = decode_state_dict(PBS_JOB_STATES, job_state)
			elif (job_state == "F"):
				exit_code = job_data["Jobs"][job_id]["Exit_status"]
				
				if (exit_code in PBS_EXIT_CODES):
					result = decode_state_dict(PBS_EXIT_CODES, exit_code)
				elif (1 <= exit_code < 128):
					result = decode_state_dict(PBS_EXIT_CODES, "failed", exit_code=exit_code)
				elif (exit_code >= 128):
					result = decode_state_dict(PBS_EXIT_CODES, "killed", signal=exit_code % 128)
				else:
					result = decode_state_dict(PBS_EXIT_CODES, "default", exit_code=exit_code)
			
			else:
				result = decode_state_dict(PBS_JOB_STATES, "default", job_status=job_state)
		else:
			result = decode_state_dict(PBS_JOB_STATES, "unknown", job_id=job_id)
	
	except Exception as ex:
		result = decode_state_dict(PBS_JOB_STATES, "error", job_id=job_id, ex=ex)
	
	return result

if __name__ == "__main__":
	# Parse arguments
	parser = argparse.ArgumentParser(description='Query the cluster job status for snakemake')
	parser.add_argument('job_id', type=str)
	parser.add_argument('--log-stderr', type=str,
		choices=[name for name, member in JobState.__members__.items()],
		default='success',
		help='Print a message to stderr each time a job has this or a ' +
		'more severe state. E.g. the default value of "success" means: ' +
		'Log all queries for jobs that have state "success" or something ' +
		'more severe, e.g. "failed", but ignore jobs that are still "running".')
	parser.add_argument('--log-file', type=str,
		choices=[name for name, member in JobState.__members__.items()],
		default='success',
		help='Print a message to --log-file-name each time a job has this or a ' +
		'more severe state. E.g. the default value of "success" means: ' +
		'Log all queries for jobs that have state "success" or something ' +
		'more severe, e.g. "failed", but ignore jobs that are still "running".')
	parser.add_argument('--log-file-name', type=str, default='clusterStatus.log',
		help='Specify the filename where messages for states specified in ' +
		'"--log-file" should be stored. Default: clusterStatus.log')
	
	args = parser.parse_args()
	
	# Map job states to log levels
	jobStateLogLevels = {
		'running': logging.DEBUG,
		'success': logging.INFO,
		'failed': logging.WARNING
	}
	
	# Setup logger
	streamHandler = StreamHandler(stream=sys.stderr)
	streamHandler.setLevel(jobStateLogLevels[args.log_stderr])
	
	fileHandler = FileHandler(args.log_file_name)
	fileHandler.setLevel(jobStateLogLevels[args.log_file])
	
	logging.basicConfig(
		level = logging.DEBUG,
		format = '[%(asctime)-15s] %(message)s',
		handlers = [streamHandler, fileHandler]
	)
	logger = logging.getLogger()
	
	job_data = query_qstat(args.job_id)
	(state, exit_code_name, exit_code_description) = decode_job_status(args.job_id, job_data)
	
	logger.log(jobStateLogLevels[state.name], "State for {job_id}: {state}. Exit code: {exit_code_name} ({exit_code_description})".format(
		job_id = args.job_id,
		state = state.name,
		exit_code_name = exit_code_name,
		exit_code_description = exit_code_description
	))
	
	print(state.name)

#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import re

from snakemake.utils import read_job_properties

defaultArgs = {
	"modules": ["Snakemake"]
}

parser=argparse.ArgumentParser(prefix_chars='+')
parser.add_argument("++depend", help="Space separated list of ids for jobs this job should depend on.")
parser.add_argument("++verbose", help="Print qsub command line to stderr before job submission.", action="store_true")
parser.add_argument("++dry-run", help="Don't execute anything. Debug purposes only.", action="store_true")
parser.add_argument("qsub_args", nargs="*")
parser.add_argument("jobscript")

args = parser.parse_args()

if (args.verbose):
	print("pbs-submit.py received the following args:\n", args, file=sys.stderr)

try:
	job_properties = read_job_properties(args.jobscript)
except Exception as e:
	print("Could not read job properties:", e, file=sys.stderr)

try:
	defaultArgs["modules"].extend(job_properties["cluster"]["modules"])
except Exception as e:
	print("Could not merge clusterArgs:", e, file=sys.stderr)

module_string = ""
for module in defaultArgs["modules"]:
	module_string += "module load {}\n".format(module)

try:
	with open(args.jobscript, "r") as f:
		jobscript_content = f.read()

	jobscript_content = re.sub('# <modules>', module_string, jobscript_content)
	
	if (args.verbose):
		print("Modified jobscript:\n", jobscript_content, file=sys.stderr)

	if not args.dry_run:
		with open(args.jobscript, "w") as f:
			f.write(jobscript_content)
except Exception as e:
	print("Could not read or modify jobscript:", e, file=sys.stderr)

depend = ""

if args.depend:
	for m in args.depend.split(" "):
		depend = depend + ":" + m
	depend = " -W \"depend=afterok" + depend + "\""

cmd = "qsub {} {} {}".format(depend, " ".join(args.qsub_args), args.jobscript)

if (args.verbose):
	print("Will submit this jobscript to the cluster using the following qsub command line:", file=sys.stderr)
	print(cmd, file=sys.stderr)

if not args.dry_run:
	try:
		res = subprocess.run(cmd, check=True, shell=True, stdout=subprocess.PIPE)
	except subprocess.CalledProcessError as e:
		raise e

res = res.stdout.decode()
print(res.strip())

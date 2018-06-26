import os
import sys

condor_dir = "/afs/cern.ch/user/d/dryu/workspace/private/DQM/Studies/SiPMGain/CMSSW_10_2_0_pre3_SiPMgain/src/HCALPFG/SiPMGains/test/condor"

def submit_single_run(run):
	command = "cmsRun sipm_gain_cfg.py run={} threads=8".format(run)
	with open("{}/run_{}.sh".format(condor_dir, run), 'w') as f:
		f.write("#!/bin/bash\n")
		f.write(command + "\n")



if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Submit many SiPM gain jobs')
	parser.add_argument('--run_file', type=str, help="Txt file with list of runs")
	parser.add_argument('--job_name', type=str, help="Job name")
	parser.add_argument('--events_per_run', type=int, default=-1, help="Number of events per run to process")
	args = parser.parse_args()

	runs = []
	with open(args.run_file, 'r') as f:
		for line in f:
			runs.append(line.strip())
	print "Run list: ",
	print runs

	with open("{}/run_{}.sh".format(condor_dir, args.job_name), 'w') as run_script:
		run_script.write("#!/bin/bash\n")
		run_script.write("#voms-proxy-init -voms cms\n")
		run_script.write("runs=( " + " ".join(runs) + " )\n")
		run_script.write("cmsRun $CMSSW_BASE/src/HCALPFG/SiPMGains/test/sipm_gain_cfg.py threads=8 processEvents={} run=${{runs[$1]}}\n".format(args.events_per_run))

	with open("{}/csub_{}.sh".format(condor_dir, args.job_name), 'w') as csub_script:
		csub_script.write("#!/bin/bash\n")
		csub_script.write("csub {}/run_{}.sh --cmssw --no_retar -n {}\n".format(condor_dir, args.job_name, len(runs)))

	cwd = os.getcwd()
	os.chdir(condor_dir)
	os.system("source {}/csub_{}.sh".format(condor_dir, args.job_name))
	os.chdir(cwd)
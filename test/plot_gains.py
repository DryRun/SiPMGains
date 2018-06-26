import os
import sys
from MyTools.RootUtils.seaborn_colors import SeabornColors
import MyTools.RootUtils.canvas_helpers as canvas_helpers
from ROOT import *
gROOT.SetBatch(True)
canvas_helpers.set_canvas_style()
seaborn_colors = SeabornColors()
seaborn_colors.load_palette("Reds_d")

def channel2string(channel):
	return "{}_ieta{}_iphi{}_depth{}".format(channel[0], channel[1], channel[2], channel[3])

if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser(description='Submit many SiPM gain jobs')
	parser.add_argument('--run_file', type=str, help="Path to txt file with list of runs")
	parser.add_argument('--run_range', type=int, nargs=2, help="Min and max runs")
	parser.add_argument('--tag', type=str, help="Tag for plot saving")
	args = parser.parse_args()

	# Load run lists
	runs = []
	with open(args.run_file, 'r') as f:
		for line in f:
			run = int(line.strip())
			if args.run_range:
				if args.run_range[0] <= run and run <= args.run_range[1]:
					runs.append(run)
			else:
				runs.append(run)
	print "Run list: ",
	print runs

	# Load data
	gains = {}
	gains_kernel = {}
	dgains = {}
	dgains_kernel = {}
	hash2channel = {}
	channel_runs = {}
	for run in runs:
		f = TFile("~/workspace/private/DQM/Studies/SiPMGain/CMSSW_10_2_0_pre3_SiPMgain/src/HCALPFG/SiPMGains/test/results/sipmgains_{}.root".format(run), "READ")
		if not f.IsOpen():
			continue
		t = f.Get("sipmGainAnalysis/gains")
		if not t:
			continue
		for i in xrange(t.GetEntries()):
			t.GetEntry(i)
			channel = (t.subdet.Data(), t.ieta, t.iphi, t.depth)
			channel_hash = hash(channel)
			if not channel_hash in hash2channel:
				hash2channel[channel_hash] = channel
				channel_runs[channel_hash] = []
			channel_runs[channel_hash].append(run)
			if not channel_hash in gains:
				gains[channel_hash] = {}
				dgains[channel_hash] = {}
				gains_kernel[channel_hash] = {}
				dgains_kernel[channel_hash] = {}
			gains[channel_hash][run] = t.gain
			dgains[channel_hash][run] = t.dgain
			gains_kernel[channel_hash][run] = t.gain_kernel
			dgains_kernel[channel_hash][run] = t.dgain_kernel

	# Convert to graphs
	h_gains = {}
	h_gains_kernel = {}
	h_gains_good = {}
	h_gains_kernel_good = {}
	fit_gains = {}
	fit_gains_kernel = {}
	for channel_hash in hash2channel.keys():
		print channel2string(hash2channel[channel_hash])
		h_gains[channel_hash] = TH1D("h_gains_{}".format(channel_hash), "Gains", len(channel_runs[channel_hash]), 0., len(channel_runs[channel_hash]))
		h_gains_kernel[channel_hash] = TH1D("h_gains_kernel_{}".format(channel_hash), "Gains", len(channel_runs[channel_hash]), 0., len(channel_runs[channel_hash]))
		h_gains_good[channel_hash] = TH1D("h_gains_good_{}".format(channel_hash), "Gains", len(channel_runs[channel_hash]), 0., len(channel_runs[channel_hash]))
		h_gains_kernel_good[channel_hash] = TH1D("h_gains_kernel_good_{}".format(channel_hash), "Gains", len(channel_runs[channel_hash]), 0., len(channel_runs[channel_hash]))
		channel_runs[channel_hash].sort()
		for i, run in enumerate(channel_runs[channel_hash]):
			bin = i + 1
			h_gains[channel_hash].SetBinContent(bin, gains[channel_hash][run])
			h_gains[channel_hash].SetBinError(bin, dgains[channel_hash][run])
			h_gains_kernel[channel_hash].SetBinContent(bin, gains_kernel[channel_hash][run])
			h_gains_kernel[channel_hash].SetBinError(bin, dgains_kernel[channel_hash][run])

			if (43. - 10. < gains[channel_hash][run]) and (gains[channel_hash][run] < 43. + 10.):
				h_gains_good[channel_hash].SetBinContent(bin, gains[channel_hash][run])
				h_gains_good[channel_hash].SetBinError(bin, dgains[channel_hash][run])
			if (43. - 10. < gains_kernel[channel_hash][run]) and (gains_kernel[channel_hash][run] < 43. + 10.):
				h_gains_kernel_good[channel_hash].SetBinContent(bin, gains_kernel[channel_hash][run])
				h_gains_kernel_good[channel_hash].SetBinError(bin, dgains_kernel[channel_hash][run])

			h_gains[channel_hash].GetXaxis().SetBinLabel(bin, str(run))
			h_gains_kernel[channel_hash].GetXaxis().SetBinLabel(bin, str(run))
			h_gains_good[channel_hash].GetXaxis().SetBinLabel(bin, str(run))
			h_gains_kernel_good[channel_hash].GetXaxis().SetBinLabel(bin, str(run))

		# Fit good points with line
		fit_gains[channel_hash] = TF1("fit_gains_{}".format(channel_hash), "[0]+[1]*x", h_gains_good[channel_hash].GetXaxis().GetXmin(), h_gains_good[channel_hash].GetXaxis().GetXmax())
		fit_gains[channel_hash].SetParameter(0, 43.)
		fit_gains[channel_hash].SetParLimits(0, 43. - 15., 43. + 15.)
		fit_gains[channel_hash].SetParameter(1, 0.)
		fit_gains[channel_hash].SetParLimits(1, -1., 1.)
		h_gains_good[channel_hash].Fit(fit_gains[channel_hash], "QWR0")

		fit_gains_kernel[channel_hash] = TF1("fit_gains_kernel_{}".format(channel_hash), "[0]+[1]*x", h_gains_kernel_good[channel_hash].GetXaxis().GetXmin(), h_gains_kernel_good[channel_hash].GetXaxis().GetXmax())
		fit_gains_kernel[channel_hash].SetParameter(0, 43.)
		fit_gains_kernel[channel_hash].SetParLimits(0, 43. - 15., 43. + 15.)
		fit_gains_kernel[channel_hash].SetParameter(1, 0.)
		fit_gains_kernel[channel_hash].SetParLimits(1, -1., 1.)
		h_gains_kernel_good[channel_hash].Fit(fit_gains_kernel[channel_hash], "QWR0")

		c_gain = TCanvas("c_gain_{}_{}".format(channel2string(hash2channel[channel_hash]), args.tag), "Gain", 1200, 600)
		c_gain.SetTicks()
		h_gains[channel_hash].SetMarkerStyle(20)
		h_gains[channel_hash].SetMarkerSize(1)
		h_gains[channel_hash].SetMarkerColor(kRed)
		h_gains[channel_hash].SetLineColor(kRed)
		h_gains[channel_hash].GetXaxis().SetTitle("Run")
		h_gains[channel_hash].GetYaxis().SetTitle("Gain [fC/pe]")
		h_gains[channel_hash].SetMinimum(0.)
		h_gains[channel_hash].SetMaximum(60.)
		h_gains[channel_hash].Draw("p")
		fit_gains[channel_hash].SetLineColor(seaborn_colors.get_root_color("Reds_d", 2))
		fit_gains[channel_hash].Draw("same")
		h_gains_good[channel_hash].SetMarkerStyle(20)
		h_gains_good[channel_hash].SetMarkerSize(1)
		h_gains_good[channel_hash].SetMarkerColor(kBlack)
		h_gains_good[channel_hash].Draw("p same")
		canvas_helpers.draw_text("Fit = ({:.2f}#pm{:.2f}) + ({:.2f}#pm{:.2f})*x".format(fit_gains[channel_hash].GetParameter(0), fit_gains[channel_hash].GetParError(0), fit_gains[channel_hash].GetParameter(1), fit_gains[channel_hash].GetParError(1)), 
			x=0.15, y=0.8, size_modifier=0.5)
		c_gain.SaveAs("~/workspace/private/DQM/Studies/SiPMGain/CMSSW_10_2_0_pre3_SiPMgain/src/HCALPFG/SiPMGains/test/results/figures/{}.pdf".format(c_gain.GetName()))

		c_gain_kernel = TCanvas("c_gain_kernel_{}_{}".format(channel2string(hash2channel[channel_hash]), args.tag), "Gain", 1200, 600)
		c_gain_kernel.SetTicks()
		h_gains_kernel[channel_hash].SetMarkerStyle(20)
		h_gains_kernel[channel_hash].SetMarkerSize(1)
		h_gains_kernel[channel_hash].SetMarkerColor(kRed)
		h_gains_kernel[channel_hash].SetLineColor(kRed)
		h_gains_kernel[channel_hash].GetXaxis().SetTitle("Run")
		h_gains_kernel[channel_hash].GetYaxis().SetTitle("Gain [fC/pe]")
		h_gains_kernel[channel_hash].SetMinimum(0.)
		h_gains_kernel[channel_hash].SetMaximum(60.)
		h_gains_kernel[channel_hash].Draw("p")
		fit_gains_kernel[channel_hash].SetLineColor(seaborn_colors.get_root_color("Reds_d", 2))
		fit_gains_kernel[channel_hash].Draw("same")
		h_gains_kernel_good[channel_hash].SetMarkerStyle(20)
		h_gains_kernel_good[channel_hash].SetMarkerSize(1)
		h_gains_kernel_good[channel_hash].SetMarkerColor(kBlack)
		h_gains_kernel_good[channel_hash].Draw("p same")
		canvas_helpers.draw_text("Fit = ({:.2f}#pm{:.2f}) + ({:.2f}#pm{:.2f})*x".format(fit_gains_kernel[channel_hash].GetParameter(0), fit_gains_kernel[channel_hash].GetParError(0), fit_gains_kernel[channel_hash].GetParameter(1), fit_gains_kernel[channel_hash].GetParError(1)), 
			x=0.15, y=0.8, size_modifier=0.5)
		c_gain_kernel.SaveAs("~/workspace/private/DQM/Studies/SiPMGain/CMSSW_10_2_0_pre3_SiPMgain/src/HCALPFG/SiPMGains/test/results/figures/{}.pdf".format(c_gain_kernel.GetName()))		
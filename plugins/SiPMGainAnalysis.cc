#include "HCALPFG/HcalTupleMaker/interface/SiPMGainAnalysis.h"

#include "FWCore/Framework/interface/ConstProductRegistry.h"
#include "FWCore/Framework/interface/ProductSelector.h"
#include "FWCore/Framework/interface/ProductSelectorRules.h"
#include "DataFormats/Provenance/interface/SelectedProducts.h"
#include "Math/LorentzVector.h"
#include "Math/Vector3D.h"

#include <map>
#include "boost/foreach.hpp"
#include <TBranch.h>
#include <TLorentzVector.h>

// --------------------------------------------------------------------------------------------------------
// Updated using the example here:
//		https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/CalibTracker/SiStripCommon/plugins/ShallowTree.cc
//		https://github.com/cms-sw/cmssw/blob/CMSSW_8_1_X/CalibTracker/SiStripCommon/interface/ShallowTree.h
// --------------------------------------------------------------------------------------------------------

SiPMGainAnalysis::SiPMGainAnalysis(const edm::ParameterSet& ps) {
	_tagHE = ps.getUntrackedParameter<edm::InputTag>("tagHE", edm::InputTag("hcalDigis"));
	_tokHE = consumes<QIE11DigiCollection>(_tagHE);	
	_subdet_to_string[HcalBarrel]  = "HB";
	_subdet_to_string[HcalEndcap]  = "HE";
	_subdet_to_string[HcalForward] = "HF";
	_subdet_to_string[HcalOuter]   = "HO";
}

SiPMGainAnalysis::~SiPMGainAnalysis() {
}

std::unique_ptr<SiPMGainAnalysis::ChannelHistogramMap> SiPMGainAnalysis::beginStream(edm::StreamID sid) const {
	return std::unique_ptr<SiPMGainAnalysis::ChannelHistogramMap>(new SiPMGainAnalysis::ChannelHistogramMap());
}

std::shared_ptr<SiPMGainAnalysis::ChannelHistogramMap> SiPMGainAnalysis::globalBeginRunSummary(edm::Run const& run, edm::EventSetup const& es) const {
	std::shared_ptr<SiPMGainAnalysis::ChannelHistogramMap> returnValue(new SiPMGainAnalysis::ChannelHistogramMap());
	edm::ESHandle<HcalDbService> dbs;
	es.get<HcalDbRecord>().get(dbs);
	const HcalElectronicsMap * emap = dbs->getHcalMapping();
	for (auto& it_gdid : emap->allPrecisionId()) {
		if (!it_gdid.isHcalDetId()) {
			continue;
		}
		HcalDetId did(it_gdid.rawId());
		if (did.subdet() == HcalEndcap) {
			char hname[50];
			sprintf(hname, "h_fc_HE_ieta%d_iphi%d_depth%d_run%d", did.ieta(), did.iphi(), did.depth(), run.runAuxiliary().run());
			returnValue->insert(std::pair<uint32_t, TH1F*>(did.rawId(), new TH1F(hname, hname, 500, 0., 500.)));
		}
	}
	return returnValue;
}

void SiPMGainAnalysis::streamBeginRun(edm::StreamID sid, edm::Run const& run, edm::EventSetup const& es) const {
	// Book stream histograms
	edm::ESHandle<HcalDbService> dbs;
	es.get<HcalDbRecord>().get(dbs);
	const HcalElectronicsMap * emap = dbs->getHcalMapping();
	for (auto& it_gdid : emap->allPrecisionId()) {
		if (!it_gdid.isHcalDetId()) {
			continue;
		}
		HcalDetId did(it_gdid.rawId());
		if (did.subdet() == HcalEndcap) {
			char hname[50];
			sprintf(hname, "h_fc_HE_ieta%d_iphi%d_depth%d_run%d", did.ieta(), did.iphi(), did.depth(), run.runAuxiliary().run());
			streamCache(sid)->insert(std::pair<uint32_t, TH1F*>(did.rawId(), new TH1F(hname, hname, 500, 0., 500.)));
		}
	}
}


void SiPMGainAnalysis::beginJob() {
 	_tree = _fs->make<TTree>("gains", "gains");
 	_tree->Branch("run", &(_gain_container.run));
 	_tree->Branch("subdet", &(_gain_container.subdet));
 	_tree->Branch("ieta", &(_gain_container.ieta));
 	_tree->Branch("iphi", &(_gain_container.iphi));
 	_tree->Branch("depth", &(_gain_container.depth));
 	_tree->Branch("gain", &(_gain_container.gain));
 	_tree->Branch("dgain", &(_gain_container.dgain));
}

void SiPMGainAnalysis::analyze(edm::StreamID sid, const edm::Event& event, const edm::EventSetup& es) const {
	// Load digis
	edm::Handle<QIE11DigiCollection> che_qie11;
	if (!event.getByToken(_tokHE, che_qie11)) {
		std::cout << "Collection QIE11DigiCollection isn't available" << _tagHE.label() << " " << _tagHE.instance() << std::endl;
		exit(1);
	}

	edm::ESHandle<HcalDbService> dbs;
	es.get<HcalDbRecord>().get(dbs);

	for (QIE11DigiCollection::const_iterator it=che_qie11->begin(); it!=che_qie11->end();
		++it)
	{
		const QIE11DataFrame digi = static_cast<const QIE11DataFrame>(*it);
		HcalDetId const& did = digi.detid();
		if (did.subdet() != HcalEndcap) {
			continue;
		}
		CaloSamples digi_fC = hcaldqm::utilities::loadADC2fCDB<QIE11DataFrame>(dbs, did, digi);
		for (int i=0; i<digi.samples(); i++) {
			double q = hcaldqm::utilities::adc2fCDBMinusPedestal<QIE11DataFrame>(dbs, digi_fC, did, digi, i);
			//std::cout << "[debug] did " << did << ", i=" << i << ", q=" << q << std::endl;
			//_spe_histograms[did.rawId()].Fill(q);
			if (streamCache(sid)->find(did.rawId()) != streamCache(sid)->end()) {
				streamCache(sid)->at(did.rawId())->Fill(q);
			} else {
				//std::cout << "[debug] WARNING : Did " << did << " is not in the streamCache object!" << std::endl;
			}
		}
	}
}

void SiPMGainAnalysis::streamEndRunSummary(edm::StreamID sid, edm::Run const& run, edm::EventSetup const& es, SiPMGainAnalysis::ChannelHistogramMap* totalGainHists) const {
	//Add the values seen in this Stream to the total for this Run
	for (auto& it_stream : (*streamCache(sid))) {
		SiPMGainAnalysis::ChannelHistogramMap::iterator it_global = totalGainHists->find(it_stream.first);
		if (it_global != totalGainHists->end()) {
			it_global->second->Add(it_stream.second);
		}
	}
	//now clear in order to be ready for the next Run
	streamCache(sid)->clear();
}

void SiPMGainAnalysis::globalEndRunSummary(edm::Run const& run, edm::EventSetup const& es, SiPMGainAnalysis::ChannelHistogramMap* totalGainHists) const {
	// Test: just print the map
	std::cout << "[SiPMGainAnalysis::globalEndRunSummary] DEBUG : Printing results" << std::endl;
	int counter = 0;
	for (auto& it : (*totalGainHists)) {
		std::cout << "[SiPMGainAnalysis::globalEndRunSummary] DEBUG : \t" << HcalDetId(it.first) << " => " << it.second << " , integral " << it.second->Integral() << std::endl; 
		it.second->Print("all");
		++counter;
		if (counter > 10) break;
	}

	// Do gain fits
	for (auto& it : (*totalGainHists)) {
		HcalDetId did(it.first);
		TH1F* hist = it.second;
		//if (hist->Integral() > 0) {
		//	_fs->cd();
		//	hist->Write();
		//}
		std::vector<int> peak_bins;
		TSpectrum* spectrum = new TSpectrum(15);
		int npeaks = spectrum->Search(hist, 1., "", 0.0015);
		if (npeaks <= 1) {
			//std::cerr << "[SiPMGainAnalysis::globalEndRunSummary] WARNING : npeaks=" << npeaks << " for channel " << did << ". hist integral = " << hist->Integral() << std::endl;
			continue;
		}
		double *peak_x = spectrum->GetPositionX();
		for (int i_peak = 0; i_peak < npeaks; ++i_peak) {
			double this_x = peak_x[i_peak];
			peak_bins.push_back(hist->GetXaxis()->FindBin(this_x));
		}
		std::sort(peak_bins.begin(), peak_bins.end());
		double gain_est = hist->GetBinCenter(peak_bins[1]) - hist->GetBinCenter(peak_bins[0]);

		std::map<int, TF1*> peak_fits;
		for (int i_peak = 0; i_peak <= 1; ++i_peak) {
			peak_fits[i_peak] = new TF1("peakfit" + i_peak, "gaus", hist->GetBinCenter(peak_bins[i_peak]) - gain_est / 3., hist->GetBinCenter(peak_bins[i_peak] + gain_est / 3.));
			hist->Fit(peak_fits[i_peak], "QR0");
		}

		_gain_container.run    = run.runAuxiliary().run();
		_gain_container.subdet = _subdet_to_string.at(did.subdet());
		_gain_container.ieta   = did.ieta();
		_gain_container.iphi   = did.iphi();
		_gain_container.depth  = did.depth();
		_gain_container.gain   = peak_fits[1]->GetParameter(1) - peak_fits[0]->GetParameter(0);
		_gain_container.dgain  = TMath::Sqrt(TMath::Power(peak_fits[0]->GetParError(1), 2) + TMath::Power(peak_fits[1]->GetParError(1), 2));

		_tree->Fill();

		delete spectrum;
		peak_fits.clear();
	}
}


void SiPMGainAnalysis::endRun(edm::Run const& run, edm::EventSetup const& es) {
/*
	// Do fits
	for (auto& it : _spe_histograms) {
		HcalDetId did(it.first);
		TH1F* hist = &(it.second);
		std::vector<int> peak_bins;
		spectrum = TSpectrum(15);
		int npeaks = spectrum.Search(hist, 1., "", 0.0015);
		double *peak_x = spectrum.GetPositionX();
		for (int i_peak = 0; i_peak < npeaks; ++i_peak) {
			double this_x = peak_x[i_peak];
			peak_bins.append(hist->GetXaxis().FindBin(this_x));
		}
		peak_bins.sort();
		double gain_est = hist.GetBinCenter(peak_bins[1]) - hist.GetBinCenter(peak_bins[0]);

		std::map<int, TF1*> peak_fits;
		for (int i_peak = 0; i_peak <= 1; ++i_peak) {
			peak_fits[i_peak] = new TF1("peakfit" + i_peak, "gaus", hist->GetBinCenter(peak_bins[i_peak]) - gain_est / 3., hist->GetBinCenter(peak_bins[i_peak] + gain_est / 3.));
			hist->Fit(peak_fits[i_peak], "QR0");
		}

		double gain = peak_fits[1]->GetParameter(1) - peak_fits[0]->GetParameter(0);
		double dgain = TMath::Sqrt(TMath::Power(peak_fits[0]->GetParError(1), 2) + TMath::Power(peak_fits[1]->GetParError(1), 2));

		// Fill TTree
		containers_int["run"] = run;
		containers_tstring["subdet"] = subdet_to_string[did.subdet()];

	}
*/
}

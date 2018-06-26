#ifndef SiPMGainAnalysis_h
#define SiPMGainAnalysis_h

#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "EventFilter/HcalRawToDigi/interface/HcalUnpacker.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalGenericDetId.h"
#include "CondFormats/HcalObjects/interface/HcalElectronicsMap.h"

#include "DQM/HcalCommon/interface/Utilities.h"

#include <string>
#include <map>
#include <vector>
#include <TTree.h>
#include <TSpectrum.h>
#include <TF1.h>
#include <TMath.h>

typedef std::map<uint32_t, TH1F*> ChannelHistogramMap;

class SiPMGainAnalysisStreamData {
public:
	ChannelHistogramMap _hist_fc; // Filled with raw SumQ values
	ChannelHistogramMap _hist_fc_kernel; // Filled with SumQ gaussian kernels for ADC granularity
	ChannelHistogramMap _hist_adc;
	TF1* _fgaus;
	int _nevents;
	inline void add(const SiPMGainAnalysisStreamData* data) {
		for (auto& it : _hist_fc) {
			if (data->_hist_fc.find(it.first) != data->_hist_fc.end()) {
				it.second->Add(data->_hist_fc.at(it.first));
			}
		}
		for (auto& it : _hist_fc_kernel) {
			if (data->_hist_fc_kernel.find(it.first) != data->_hist_fc_kernel.end()) {
				it.second->Add(data->_hist_fc_kernel.at(it.first));
			}
		}
		for (auto& it : _hist_adc) {
			if (data->_hist_adc.find(it.first) != data->_hist_adc.end()) {
				it.second->Add(data->_hist_adc.at(it.first));
			}
		}
		_nevents += data->_nevents;
	}

	inline void reset() {
		this->_hist_fc.clear();
		this->_hist_fc_kernel.clear();
		this->_hist_adc.clear();
		if (_fgaus) {
			delete _fgaus;
			_fgaus = 0;
		}
	}
};

class SiPMGainAnalysis :  
	public edm::global::EDAnalyzer<edm::StreamCache<SiPMGainAnalysisStreamData>, edm::RunSummaryCache<SiPMGainAnalysisStreamData>> {
public:
	edm::InputTag		_tagHE;
	edm::EDGetTokenT<QIE11DigiCollection> _tokHE;

public:
	explicit SiPMGainAnalysis(const edm::ParameterSet& ps);// : pset(iConfig) {}
	~SiPMGainAnalysis();

 private:
	
	void beginJob();
	void beginRun(edm::Run const&, edm::EventSetup const&) {}
	std::shared_ptr<SiPMGainAnalysisStreamData> globalBeginRunSummary(edm::Run const&, edm::EventSetup const&) const;
	std::unique_ptr<SiPMGainAnalysisStreamData> beginStream(edm::StreamID) const;
	void streamBeginRun(edm::StreamID, edm::Run const&, edm::EventSetup const&) const;
	void analyze(edm::StreamID, const edm::Event&, const edm::EventSetup&) const;
	void streamEndRunSummary(edm::StreamID, edm::Run const&, edm::EventSetup const&, SiPMGainAnalysisStreamData*) const;
	void globalEndRunSummary(edm::Run const&, edm::EventSetup const&, SiPMGainAnalysisStreamData*) const;

	void endJob(){}
	void endRun(edm::Run const&, edm::EventSetup const&){}
	
	std::map<TString, double> fitGain(TH1F* hist, TF1** peak_fit, bool debug=false) const;

	edm::Service<TFileService> _fs;
	//std::map<uint32_t, TH1F*> _spe_histograms;
	std::map<HcalSubdetector, TString> _subdet_to_string;

	struct GainContainer {
		int run                = 0;
		int nevents            = 0;
		TString subdet         = "";
		int ieta               = 0;
		int iphi               = 0;
		int depth              = 0;
		double gain            = 0.;
		double gain_est        = 0.;
		double dgain           = 0.;
		double mean1           = 0.;
		double mean2           = 0.;
		double gain_kernel     = 0.;
		double gain_est_kernel = 0.;
		double dgain_kernel    = 0.;
		double mean1_kernel    = 0.;
		double mean2_kernel    = 0.;
	};

	mutable GainContainer _gain_container;
	TTree *_tree;


protected:
	edm::ESHandle<HcalDbService> _dbService;
	
 public:
	static const double adc2fC[];
};

#endif

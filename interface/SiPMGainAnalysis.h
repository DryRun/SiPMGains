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

class SiPMGainAnalysis :  
	public edm::global::EDAnalyzer<edm::StreamCache<std::map<uint32_t, TH1F*>>, edm::RunSummaryCache<std::map<uint32_t, TH1F*>>> {
public:
	typedef std::map<uint32_t, TH1F*> ChannelHistogramMap;
	edm::InputTag		_tagHE;
	edm::EDGetTokenT<QIE11DigiCollection> _tokHE;
 private:
	
	void beginJob();
	void beginRun(edm::Run const&, edm::EventSetup const&) {}
	std::shared_ptr<ChannelHistogramMap> globalBeginRunSummary(edm::Run const&, edm::EventSetup const&) const;
	std::unique_ptr<ChannelHistogramMap> beginStream(edm::StreamID) const;
	void streamBeginRun(edm::StreamID, edm::Run const&, edm::EventSetup const&) const;
	void analyze(edm::StreamID, const edm::Event&, const edm::EventSetup&) const;
	void streamEndRunSummary(edm::StreamID, edm::Run const&, edm::EventSetup const&, ChannelHistogramMap*) const;
	void globalEndRunSummary(edm::Run const&, edm::EventSetup const&, ChannelHistogramMap*) const;

	void endJob(){}
	void endRun(edm::Run const&, edm::EventSetup const&);
	
	edm::Service<TFileService> _fs;
	//std::map<uint32_t, TH1F*> _spe_histograms;
	std::map<HcalSubdetector, TString> _subdet_to_string;

	struct GainContainer {
		int run = 0;
		TString subdet = "";
		int ieta = 0;
		int iphi = 0;
		int depth = 0;
		double gain = 0.;
		double dgain = 0.;
	};

	mutable GainContainer _gain_container;
	TTree *_tree;
	
protected:
	edm::ESHandle<HcalDbService> _dbService;
	
 public:
	explicit SiPMGainAnalysis(const edm::ParameterSet& ps);// : pset(iConfig) {}
	~SiPMGainAnalysis();
	
};

#endif

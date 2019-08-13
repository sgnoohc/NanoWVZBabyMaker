#ifndef analysis_h
#define analysis_h

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "rooutil.h"
#include "Nano.h"

class AnalysisConfig {

public:

    // TString that holds the input file list (comma separated)
    TString input_file_list_tstring;

    // TString that holds the name of the TTree to open for each input files
    TString input_tree_name;

    // Output TFile
    TFile* output_tfile;

    // Output TTree
    TTree* output_ttree;

    // Number of events to loop over
    int n_events;

    // Jobs to split (if this number is positive, then we will skip certain number of events)
    // If there are N events, and was asked to split 2 ways, then depending on job_index, it will run over first half or latter half
    int nsplit_jobs;

    // Job index (assuming nsplit_jobs is set, the job_index determine where to loop over)
    int job_index;

    // Debug boolean
    bool debug;

    // Debug boolean
    bool is_data;

    // TChain that holds the input TTree's
    TChain* events_tchain;

    // Custom Looper object to facilitate looping over many files
    RooUtil::Looper<Nano> looper;

    // Custom Cutflow framework
    RooUtil::Cutflow cutflow;

    // Custom Histograms object compatible with RooUtil::Cutflow framework
    RooUtil::Histograms histograms;

    // Custom TTreeX object to facilitate output ttree
    RooUtil::TTreeX* tx;

    // Custom processor object that will run per event to process events in a modular way
    RooUtil::Processor* processor;

};

extern AnalysisConfig ana;

#endif

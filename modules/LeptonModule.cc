#include "LeptonModule.h"

//==============================================================================================================
//
// Lepton Module
//
//==============================================================================================================

//##############################################################################################################
void wvzModule::LeptonModule::AddOutput()
{

    tx->createBranch<vector<LV>>("lep_p4");
    tx->createBranch<vector<float>>("lep_pt");
    tx->createBranch<vector<float>>("lep_eta");
    tx->createBranch<vector<float>>("lep_phi");
    tx->createBranch<vector<float>>("lep_pfRelIso04_all");
    tx->createBranch<vector<float>>("lep_pfRelIso03_all");
    tx->createBranch<vector<float>>("lep_ip3d");
    tx->createBranch<vector<float>>("lep_sip3d");
    tx->createBranch<vector<float>>("lep_dxy");
    tx->createBranch<vector<float>>("lep_dz");
    tx->createBranch<vector<float>>("lep_deltaEtaSC");
    tx->createBranch<vector<int>>("lep_genPartFlav");
    tx->createBranch<vector<int>>("lep_idx");
    tx->createBranch<vector<int>>("lep_id");
    tx->createBranch<vector<int>>("lep_isWVZVeto");
    tx->createBranch<vector<int>>("lep_isWVZNominal");
    tx->createBranch<int>("nVlep");

}

//##############################################################################################################
void wvzModule::LeptonModule::FillOutput()
{

    // I expect either one electron or one muon so it's ok to loop over
    for (unsigned int idx = 0; idx < nt.Electron_pt().size(); ++idx)
    {
        tx->pushbackToBranch<LV>("lep_p4", nt.Electron_p4()[idx]);
        tx->pushbackToBranch<float>("lep_pt", nt.Electron_pt()[idx]);
        tx->pushbackToBranch<float>("lep_eta", nt.Electron_eta()[idx]);
        tx->pushbackToBranch<float>("lep_phi", nt.Electron_phi()[idx]);
        tx->pushbackToBranch<float>("lep_pfRelIso04_all", -999);
        tx->pushbackToBranch<float>("lep_pfRelIso03_all", nt.Electron_pfRelIso03_all()[idx]);
        tx->pushbackToBranch<float>("lep_ip3d", nt.Electron_ip3d()[idx]);
        tx->pushbackToBranch<float>("lep_sip3d", nt.Electron_sip3d()[idx]);
        tx->pushbackToBranch<float>("lep_dxy", nt.Electron_dxy()[idx]);
        tx->pushbackToBranch<float>("lep_dz", nt.Electron_dz()[idx]);
        tx->pushbackToBranch<float>("lep_deltaEtaSC", nt.Electron_deltaEtaSC()[idx]);
        tx->pushbackToBranch<int>("lep_genPartFlav", nt.Electron_genPartFlav()[idx]);
        tx->pushbackToBranch<int>("lep_idx", idx);
        tx->pushbackToBranch<int>("lep_id", nt.Electron_charge()[idx]*(-11));
        tx->pushbackToBranch<int>("lep_isWVZVeto", isPt10AnalysisVetoElectron(idx));
        tx->pushbackToBranch<int>("lep_isWVZNominal", isPt10AnalysisNominalElectron(idx));
    }

    // I expect either one electron or one muon so it's ok to loop over
    for (unsigned int idx = 0; idx < nt.Muon_pt().size(); ++idx)
    {
        tx->pushbackToBranch<LV>("lep_p4", nt.Muon_p4()[idx]);
        tx->pushbackToBranch<float>("lep_pt", nt.Muon_pt()[idx]);
        tx->pushbackToBranch<float>("lep_eta", nt.Muon_eta()[idx]);
        tx->pushbackToBranch<float>("lep_phi", nt.Muon_phi()[idx]);
        tx->pushbackToBranch<float>("lep_pfRelIso04_all", nt.Muon_pfRelIso04_all()[idx]);
        tx->pushbackToBranch<float>("lep_pfRelIso03_all", nt.Muon_pfRelIso03_all()[idx]);
        tx->pushbackToBranch<float>("lep_ip3d", nt.Muon_ip3d()[idx]);
        tx->pushbackToBranch<float>("lep_sip3d", nt.Muon_sip3d()[idx]);
        tx->pushbackToBranch<float>("lep_dxy", nt.Muon_dxy()[idx]);
        tx->pushbackToBranch<float>("lep_dz", nt.Muon_dz()[idx]);
        tx->pushbackToBranch<float>("lep_deltaEtaSC", -999);
        tx->pushbackToBranch<int>("lep_genPartFlav", nt.Muon_genPartFlav()[idx]);
        tx->pushbackToBranch<int>("lep_idx", idx);
        tx->pushbackToBranch<int>("lep_id", nt.Muon_charge()[idx]*(-13));
        tx->pushbackToBranch<int>("lep_isWVZVeto", isPt10AnalysisVetoMuon(idx));
        tx->pushbackToBranch<int>("lep_isWVZNominal", isPt10AnalysisNominalMuon(idx));
    }

    tx->sortVecBranchesByPt("lep_p4",
            {
            "lep_pt",
            "lep_eta",
            "lep_phi",
            "lep_pfRelIso04_all",
            "lep_pfRelIso03_all",
            "lep_ip3d",
            "lep_sip3d",
            "lep_dxy",
            "lep_dz",
            "lep_deltaEtaSC",
            },
            {
            "lep_idx",
            "lep_id",
            "lep_isWVZVeto",
            "lep_isWVZNominal",
            },
            {});

    tx->setBranch<int>("nVlep", passCount(tx->getBranchLazy<vector<int>>("lep_isWVZVeto")));

}

//##############################################################################################################
bool wvzModule::LeptonModule::isPt10AnalysisVetoMuon(int idx)
{
    if (!( nt.Muon_p4()[idx].pt() > 10.            )) return false;
    if (!( nt.Muon_mediumId()[idx]                 )) return false;
    if (!( fabs(nt.Muon_p4()[idx].eta()) < 2.4     )) return false;
    if (!( nt.Muon_pfRelIso04_all()[idx]  < 0.25   )) return false;
    if (fabs(nt.Muon_p4()[idx].eta()) <= 1.479)
    {
        if (!( fabs(nt.Muon_dz()[idx]) < 0.1       )) return false;
        if (!( fabs(nt.Muon_dxy()[idx]) < 0.05     )) return false;
    }
    else
    {
        if (!( fabs(nt.Muon_dz()[idx]) < 0.2       )) return false;
        if (!( fabs(nt.Muon_dxy()[idx]) < 0.1      )) return false;
    }
    if (!( fabs(nt.Muon_sip3d()[idx]) < 4. )) return false;
    return true;
}

//##############################################################################################################
bool wvzModule::LeptonModule::isPt10AnalysisVetoElectron(int idx)
{
    if (!( nt.Electron_p4()[idx].pt() > 10.            )) return false;
    if (!( nt.Electron_cutBased()[idx] >= 1            )) return false;
    if (!( fabs(nt.Electron_p4()[idx].eta()) < 2.5     )) return false;
    if (fabs(nt.Electron_deltaEtaSC()[idx] + nt.Electron_eta()[idx]) <= 1.479)
    {
        if (!( fabs(nt.Electron_dz()[idx]) < 0.1       )) return false;
        if (!( fabs(nt.Electron_dxy()[idx]) < 0.05     )) return false;
    }
    else
    {
        if (!( fabs(nt.Electron_dz()[idx]) < 0.2       )) return false;
        if (!( fabs(nt.Electron_dxy()[idx]) < 0.1      )) return false;
    }
    if (!( fabs(nt.Electron_sip3d()[idx]) < 4.         )) return false;
    return true;
}

//##############################################################################################################
bool wvzModule::LeptonModule::isPt10AnalysisNominalMuon(int idx)
{
    if (!( isPt10AnalysisVetoMuon(idx) )) return false;
    if (!( nt.Muon_pfRelIso04_all()[idx] < 0.15 )) return false;
    return true;
}

//##############################################################################################################
bool wvzModule::LeptonModule::isPt10AnalysisNominalElectron(int idx)
{
    if (!( isPt10AnalysisVetoElectron(idx) )) return false;
    if (!( nt.Electron_cutBased()[idx] >= 3 )) return false;
    return true;
}

//##############################################################################################################
int wvzModule::LeptonModule::passCount(const vector<int>& v)
{
    return std::count_if(v.begin(), v.end(), [](int i){return i > 0;});
}

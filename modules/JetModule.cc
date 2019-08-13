#include "JetModule.h"

//==============================================================================================================
//
// MET Module
//
//==============================================================================================================

void wvzModule::JetModule::AddOutput()
{

    tx->createBranch<vector<LV>>("jets_p4");
    tx->createBranch<vector<float>>("jets_btag_score");
    tx->createBranch<int>("nb");
}

void wvzModule::JetModule::FillOutput()
{

    const vector<int>& lep_id = tx->getBranchLazy<vector<int>>("lep_id");
    const vector<int>& lep_idx = tx->getBranchLazy<vector<int>>("lep_idx");
    const vector<int>& lep_isWVZVeto = tx->getBranchLazy<vector<int>>("lep_isWVZVeto");

    int nb = 0;

    for (unsigned int idx = 0; idx < nt.Jet_pt().size(); ++idx)
    {

        if (gconf.year == 2017)
        {
            // 2017 does not have ID?? (FIXME: Check this)
        }
        else
        {
            if (!( nt.Jet_jetId()[idx] & (1<<1) ))
                continue;
        }

        LV jet = nt.Jet_p4()[idx];
        float btag_score = nt.Jet_btagCSVV2()[idx];

        // Overlap remove
        for (unsigned int ilep = 0; ilep < lep_id.size(); ++ilep)
        {

            // Only overlap remove against veto leptons
            if (lep_isWVZVeto[ilep])
                continue;

            if (abs(lep_id[ilep]) == 11)
            {
                if (nt.Jet_electronIdx1()[idx] == lep_idx[ilep])
                    continue;
            }
            else
            {
                if (nt.Jet_muonIdx1()[idx] == lep_idx[ilep])
                    continue;
            }

        }

        float btag_loose_threshold = -999;
        if (gconf.year == 2016) btag_loose_threshold = 0.2217;
        if (gconf.year == 2017) btag_loose_threshold = 0.1522;
        if (gconf.year == 2018) btag_loose_threshold = 0.1241;

        // Now ask whether the jet has the kinematic requirements of a b-jet
        if (jet.pt() > 20. and fabs(jet.eta()) < 2.4)
        {
            if (btag_score >= btag_loose_threshold)
                nb++;
        }

        if (jet.pt() > 30.)
        {
            tx->pushbackToBranch<LV>("jets_p4", jet);
            tx->pushbackToBranch<float>("jets_btag_score", btag_score);
        }

    }

    tx->setBranch<int>("nb", nb);

}


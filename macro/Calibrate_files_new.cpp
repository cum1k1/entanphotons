#include <iostream>
#include <fstream>
#include "like_ivashkin_wants_it.h"

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TChain.h>
#include <TObjArray.h>
#include <TString.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TObjString.h>
#include "../calo_analysis/Readme_reader.h"
#include "ChannelEntry.h"

#define UsedScatterer 0

const int GATE_BEG = 50;
const int GATE_END = 150;
const Int_t N_zl_intervals = 2;
const Int_t Gate_end_for_main_scatterers = 71;
////////////////////////
///////////////////////Define fit constants range
Int_t l_integr_NaI = 40000;
Int_t r_integr_NaI = 300000;

///////////////////////
///////////////////////

void Calibrate_files_new()
{

	TString source_path = "/home/doc/entanglement/new_files_data/new_run/";
    #if UsedScatterer
    const Int_t total_channels = 35;
    #else
    const Int_t total_channels = 34;
    #endif

    ChannelEntry *channel_info = new ChannelEntry[total_channels];
    short_energy_ChannelEntry *short_channel_info = new short_energy_ChannelEntry[total_channels];
    //PMT_tree->SetBranchAddress(Form("channel_%d", channel_number), &channel_info[channel_number]);    
    TFile *combined_root = new TFile (source_path+"calibrated_time.root", "RECREATE");
    TTree *combined_tree = new TTree ("adc64_data","adc64_data");
    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
        (short_channel_info+channel_number)->CreateBranch(combined_tree, channel_number);
    Float_t interval_width = 1.;
	TObjArray files_names;
	TSystemDirectory dir(source_path, source_path);
	TList *files = dir.GetListOfFiles();
	if(files) 
	{
		TSystemFile *file; 
		TString fname; 
		TIter next(files); 
		while((file=(TSystemFile*)next())) 
		{           
			fname = file->GetName(); 
			if(!file->IsDirectory() && fname.EndsWith(".root") && fname != "combined.root" && fname != "calibrated.root"&& fname != "calibrated_time.root" && fname != "entangled_mini_tree.root" && fname != "decoherent_mini_tree.root") 
			{
                TChain *PMT_tree = new TChain;

				files_names.Add(new TObjString(fname));
				PMT_tree->AddFile( (source_path + fname + "/adc64_data").Data() );
                for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
                    (channel_info+channel_number)->SetBranch(PMT_tree,channel_number);

///////////////////////
///////////////////////Calculating needed values

///////////////////////
///////////////////////
                Int_t events_divider = 2;
                Int_t cut_number = events_divider;
                Int_t low_cut[total_channels] = {0}; Int_t high_cut[total_channels] = {0};
                Float_t sigma_i[total_channels] = {0}; Float_t mean_int[total_channels] = {0};
                Int_t total_entries = PMT_tree->GetEntries()/100;
                Int_t Events_for_cuts_left = PMT_tree->GetEntries()/2000;
                Int_t Events_for_cuts_right = Events_for_cuts_left + PMT_tree->GetEntries()/100;


                TH1F *peak_histo[total_channels];
                for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                {
                    TString hernya = Form("integral_in_gate_ADC_ch_%i",channel_number);
                    peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),100,l_integr_NaI,r_integr_NaI);
                }

                for (Int_t iEvent = Events_for_cuts_left; iEvent < Events_for_cuts_right; iEvent++)
                {
                    PMT_tree->GetEntry(iEvent);
                    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                    {
                        Int_t zero_level = channel_info[channel_number].Get_Zero_Level(GATE_BEG);
                        //cout<<"ch="<<channel_number<<"Amp = "<<channel_info[channel_number].Get_Amplitude(zero_level)<<" zl="<<zero_level<<endl;
                    
                        Int_t gateIntegral = channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END);
                        if (gateIntegral > 4000
                        && channel_info[channel_number].Get_Amplitude(zero_level) > 400
                        && channel_info[channel_number].Get_Amplitude(zero_level) < 60000

            #if UsedScatterer
                        && abs(channel_info[34].peak_pos-channel_info[32].peak_pos) > 10
                        && channel_info[34].amp < 400
            #endif
                        ) 
                        {
                            cout<<"ch="<<channel_number<<"Amp = "<<channel_info[channel_number].Get_Amplitude(zero_level)<<" zl="<<zero_level<<endl;

                            peak_histo[channel_number]->Fill(gateIntegral);
                        }
                    }
                }

                for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                {
                    //PMT_tree->Draw(Form("channel_%i.integral_in_gate >> %s",channel_number, hernya.Data()), Form("channel_%i.integral_in_gate > 4000",channel_number));
                    TF1 *gaus_fit = new TF1 ("gaus_func","gaus",
                    peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())-20000,
                    peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())+20000);
                    peak_histo[channel_number]->Fit("gaus_func","R");

                    TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
                    (gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
                    (gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
                    peak_histo[channel_number]->Fit("gaus_func_2","R");

                    mean_int[channel_number] = gaus_2_fit->GetParameter(1);
                    sigma_i[channel_number] = gaus_2_fit->GetParameter(2);
                    low_cut[channel_number] = gaus_2_fit->GetParameter(1)-interval_width*gaus_2_fit->GetParameter(2);
                    high_cut[channel_number] = gaus_2_fit->GetParameter(1)+interval_width*gaus_2_fit->GetParameter(2);
                    peak_histo[channel_number]->Write();
                    
                    delete peak_histo[channel_number];

                }
            ////////////////////////////////////
                TString left_name = "integral_in_gate_left_arm_scatterer";
                TString right_name = "integral_in_gate_right_arm_scatterer";
                Int_t scatterer_left_boarder_peak_search = 00000;


                peak_histo[32] = new TH1F (left_name.Data(),left_name.Data(),200,scatterer_left_boarder_peak_search,500000);
                peak_histo[33] = new TH1F (right_name.Data(),right_name.Data(),200,scatterer_left_boarder_peak_search,500000);

                    //TH1F *peak_histo[32] = nullptr;
                for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                {
                    TString hernya = Form("integral_in_gate_scatterer_ch_cut_%i",channel_number);
                    peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),200,scatterer_left_boarder_peak_search,500000);
                }


                for (Int_t iEvent = Events_for_cuts_left; iEvent < Events_for_cuts_right; iEvent++)
                {
                    PMT_tree->GetEntry(iEvent);
                    for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                    {
                        Int_t zero_level = channel_info[channel_number].Get_Zero_Level(GATE_BEG);

                        Int_t sc_number = 0;
                        if (channel_number < 16) sc_number = 32;
                        if (channel_number >= 16) sc_number = 33;
                        
                        if (channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END) > low_cut[channel_number]
                        && channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END) < high_cut[channel_number]
                        //&& abs(channel_info[channel_number].Get_peak_position(zero_level, GATE_BEG, GATE_END)-channel_info[sc_number].Get_peak_position(zero_level, GATE_BEG, GATE_END)) < 20
                        && (channel_info[channel_number].Get_peak_position(zero_level, GATE_BEG, GATE_END)-channel_info[sc_number].Get_peak_position(zero_level, GATE_BEG, GATE_END)) > 0
                        && (channel_info[channel_number].Get_peak_position(zero_level, GATE_BEG, GATE_END)-channel_info[sc_number].Get_peak_position(zero_level, GATE_BEG, GATE_END)) < 14

                        && channel_info[channel_number].Get_Amplitude(zero_level) < 60000
                        //&& channel_info[sc_number].Get_peak_position(zero_level, GATE_BEG, GATE_END) < 31000
                        

            #if UsedScatterer
                        && abs(channel_info[34].peak_pos-channel_info[32].peak_pos > 10)
                        && channel_info[34].amp < 400

            #endif
                        )
                        peak_histo[channel_number]->Fill(channel_info[sc_number].Get_Integral_in_Gate(zero_level, GATE_BEG, Gate_end_for_main_scatterers));
                    }
                }


                for (Int_t channel_number = 0; channel_number < 32; channel_number++)
                {

                    peak_histo[channel_number]->Write();
                    if (channel_number < 16) peak_histo[32]->Add(peak_histo[channel_number]);

                    if (channel_number >= 16) peak_histo[33]->Add(peak_histo[channel_number]);

                    delete peak_histo[channel_number];
                }

                for (Int_t channel_number = 32; channel_number < 34; channel_number++)
                {
                    TF1 *gaus_fit = new TF1 ("gaus_func","gaus",
                    peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())-20000,
                    peak_histo[channel_number]->GetBinCenter(peak_histo[channel_number]->GetMaximumBin())+20000);
                    peak_histo[channel_number]->Fit("gaus_func","R");

                    TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
                    (gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
                    (gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
                    peak_histo[channel_number]->Fit("gaus_func_2","R");
                    mean_int[channel_number] = gaus_2_fit->GetParameter(1);
                    peak_histo[channel_number]->Write();
                    delete peak_histo[channel_number];

                }
//////////////////////
//////////////////////Fill_clibrated_tree
                for (Long64_t entry_num = 0; entry_num < total_entries; entry_num++)
                {
                    PMT_tree->GetEntry(entry_num);
                    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
                    {
                        Int_t zero_level = channel_info[channel_number].Get_Zero_Level(GATE_BEG);
                        Int_t CH_GATE_END = GATE_END; if (channel_number == 32 || channel_number ==33) CH_GATE_END = Gate_end_for_main_scatterers;

                        Int_t energy_in_peak = 255;
                        if (channel_number == 32 || channel_number == 33) energy_in_peak = 255;
                        if (channel_number < 34) 
                        {
                            short_channel_info[channel_number].integral_in_gate = 
                            (float)energy_in_peak/mean_int[channel_number]*
                            channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, CH_GATE_END);
                            if ((entry_num + 1) % 100000 == 0) cout << channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, CH_GATE_END) << " " << short_channel_info[channel_number].integral_in_gate << endl;
                        }
                        else 
                        {if (channel_number == 34) short_channel_info[channel_number].integral_in_gate = 
                        (float)30./60000*channel_info[channel_number].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END);}
                        //Int_t inv_amp = channel_info[channel_number].Get_Amplitude(zero_level);
                        //short &wf_size = channel_info[channel_number].wf_size;
                        Short_t *wf = channel_info[channel_number].wf;

                        if ((entry_num + 1) % 100000 == 0) cout << channel_info[channel_number].Get_Amplitude(zero_level) << endl<<endl<<endl;
                        short_channel_info[channel_number].amp = channel_info[channel_number].Get_Amplitude(zero_level, GATE_BEG, CH_GATE_END);//inv_amp;

                        short_channel_info[channel_number].peak_pos = 
                        //channel_info[channel_number].Get_peak_position_gauss(zero_level, short_channel_info[channel_number].amp, GATE_BEG, CH_GATE_END);
                        channel_info[channel_number].Get_peak_position(zero_level, GATE_BEG, CH_GATE_END);
                        if (channel_info[channel_number].wf_size == 350 && (entry_num + 1) % 100000 == 0) cout << short_channel_info[channel_number].peak_pos << " " << channel_info[channel_number].Get_peak_position(zero_level, GATE_BEG, CH_GATE_END) << "; channel_number = " << channel_number<< endl;
                    }
                    combined_tree->Fill();
                    if ((entry_num + 1) % 100000 == 0) cout << (float)entry_num/total_entries*100<<"%"<<endl;
                }
                delete PMT_tree;
////////////////////////

			}
		}

	}

    
    combined_tree->Write();
    combined_root->Close();
}
        




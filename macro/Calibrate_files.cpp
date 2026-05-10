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

const int GATE_BEG = 150;
const int GATE_END = 300;


void Calibrate_files()
{

	TString source_path = "/home/daq/entanglement/tmp/";
    const Int_t total_channels = 35;

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
			if(!file->IsDirectory() && fname.EndsWith(".root") && fname != "combined.root" && fname != "calibrated.root"&& fname != "calibrated_time.root") 
			{
                TChain *PMT_tree = new TChain;

				files_names.Add(new TObjString(fname));
				PMT_tree->AddFile( (source_path + fname + "/adc64_data").Data() );
                Int_t total_entries = PMT_tree->GetEntries();
                for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
                    (channel_info+channel_number)->SetBranch(PMT_tree,channel_number);


///////////////////////
///////////////////////
    Int_t low_cut[32] = {0}; Int_t high_cut[32] = {0};
    Float_t sigma_i[32] = {0}; Float_t mean_int[34] = {0};
    Int_t Events_for_cuts = total_entries;
///
        TH1F *peak_histo[35];
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            TString hernya = Form("integral_in_gate_ADC_ch_%i",channel_number);
            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),100,100000,200000);
        }

        for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                if (channel_info[channel_number].integral_in_gate > 4000
                && channel_info[channel_number].amp > 400
                && channel_info[channel_number].amp < 65000


#if UsedScatterer
                && abs(channel_info[34].peak_pos-channel_info[32].peak_pos) > 10
                && channel_info[34].amp < 400
#endif
                ) 
                {
                    peak_histo[channel_number]->Fill(channel_info[channel_number].integral_in_gate);
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
                
                delete peak_histo[channel_number];

            }
////////////////////////////////////
    TString left_name = "integral_in_gate_left_arm_scatterer";
    TString right_name = "integral_in_gate_right_arm_scatterer";

    peak_histo[32] = new TH1F (left_name.Data(),left_name.Data(),200,1000,500000);
    peak_histo[33] = new TH1F (right_name.Data(),right_name.Data(),200,1000,500000);

        //TH1F *peak_histo[32] = nullptr;
        for (Int_t channel_number = 0; channel_number < 32; channel_number++)
        {
            TString hernya = Form("integral_in_gate_scatterer_ch_cut_%i",channel_number);
            peak_histo[channel_number] = new TH1F (hernya.Data(),hernya.Data(),200,1000,500000);
        }


        for (Int_t iEvent = 0; iEvent < Events_for_cuts; iEvent++)
        {
            PMT_tree->GetEntry(iEvent);
            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
                Int_t sc_number = 0;
                if (channel_number < 16) sc_number = 32;
                if (channel_number >= 16) sc_number = 33;
                if (channel_info[channel_number].integral_in_gate > low_cut[channel_number]
                && channel_info[channel_number].integral_in_gate < high_cut[channel_number]
                && abs(channel_info[channel_number].peak_pos-channel_info[sc_number].peak_pos) < 20
                && channel_info[channel_number].amp < 65000
                && channel_info[sc_number].amp < 65000
               

#if UsedScatterer
                && abs(channel_info[34].peak_pos-channel_info[32].peak_pos > 10)
                && channel_info[34].amp < 400

#endif
                )
                peak_histo[channel_number]->Fill(channel_info[sc_number].integral_in_gate);
            }
        }


            for (Int_t channel_number = 0; channel_number < 32; channel_number++)
            {
  
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
                for (Long64_t entry_num = 0; entry_num < PMT_tree->GetEntries(); entry_num++)
                {
                    PMT_tree->GetEntry(entry_num);
                    for (Int_t channel_number = 0; channel_number < total_channels; channel_number++)
                    {
                        Int_t energy_in_peak = 260;
                        if (channel_number == 32 || channel_number == 33) energy_in_peak = 250;
                        if (channel_number < 34) 
                        {
                            short_channel_info[channel_number].integral_in_gate = 
                            (float)energy_in_peak/mean_int[channel_number]*channel_info[channel_number].integral_in_gate;
                            //if (entry_num%100==0 && (channel_number == 33 || channel_number ==32 || channel_number == 31)) cout <<"channel_number = "<<channel_number<< "; short_in_gate = "<<short_channel_info[channel_number].integral_in_gate << "; int_in_gate = " << channel_info[channel_number].integral_in_gate<< "; mean_int = "<< mean_int[channel_number] << endl;
                        }
                        else 
                        {if (channel_number == 34) short_channel_info[channel_number].integral_in_gate = 
                        (float)30/60000*channel_info[channel_number].integral_in_gate;}
                        Int_t inv_amp = 0;
                        short &wf_size = channel_info[channel_number].wf_size;
                        Short_t *wf = channel_info[channel_number].wf; 
                        Short_t &zl = channel_info[channel_number].zero_level;
                        for (int s=0; s < wf_size; ++s) {
                            Int_t v = zl - wf[s];
                            if (v > inv_amp) inv_amp = v;
                        }
                        short_channel_info[channel_number].amp = (UShort_t) inv_amp;
                        Float_t peak_search = 0.;
                        Float_t ampl_sum = 0;
                        Int_t count_gates = wf_size;
                        if (channel_number == 33) count_gates = 240;
                        for (int s=0; s < count_gates; ++s) {
                            Int_t v = zl - wf[s];
                            if (v > inv_amp*0.1)
                            {
                                ampl_sum += (float) v;
                                peak_search+= (float) v*s;

                            }

                        }
                            peak_search /= ampl_sum;

                        //short_channel_info[channel_number].peak_pos = channel_info[channel_number].peak_pos;
                        short_channel_info[channel_number].peak_pos = (Short_t) 16.0*peak_search;
                        if (channel_info[channel_number].wf_size == 350 && (entry_num + 1) % 100000 == 0) cout << short_channel_info[channel_number].peak_pos << " " << channel_info[channel_number].peak_pos << "; channel_number = " << channel_number<< endl;
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
        




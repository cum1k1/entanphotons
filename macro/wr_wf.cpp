#include <iostream>
#include <fstream>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TString.h>
#include <stdio.h>

#include <stdlib.h>
#include <sstream>
#include "dat_reader.h"

using namespace std;


void wr_wf(TString file_path = "/entanglement_new_files/BOCHKA/Termos2/test2_azot_cs137/",    TString file_name = "TR_0_0", TString tempdir = "")
{
    const Int_t GATE_BEG = 50;
    const Int_t GATE_END = 1000;
    const Int_t total_channels = 1;
    const Int_t Event_Length = 1024;




    ecal_data channel_info[total_channels];
    for (Int_t i = 0; i < total_channels; i++) channel_info[i].Initialize();
    short_energy_ChannelEntry short_channel_info[total_channels];
    ifstream inputfile;
    inputfile.open(file_path+file_name+".txt");
    FILE *DataFile = fopen(file_path+file_name+".dat","rb");
    TFile *combined_root = new TFile (file_path+file_name+".root", "RECREATE");
    TTree *combined_tree = new TTree ("adc64_data","adc64_data");
    for (Int_t chnum = 0; chnum < total_channels; chnum++)
        (short_channel_info[chnum]).CreateBranch(combined_tree, chnum);  

    Int_t event_number=0;
    string waveform;
    Int_t coincidence_counter[4][2] = {0};
    Int_t global_counter_ev = 0;
    while(fread(channel_info[0].wf, 4, Event_Length+6, DataFile)==Event_Length+6)
    {            
            Int_t flag = 0;
            channel_info[0].wf_size = Event_Length+6;
        //stringstream ss(waveform);

            Short_t min_diff,min_time,max_diff,max_time;
            Int_t chnum = 0;

            //ss >> chnum;
            //if (chnum == 64) cout << event_number << endl;

            // stoi(waveform.substr(0,waveform.find(" ")));
            //channel_info[chnum].ReadWf(waveform.substr(waveform.find(" ")+1));
            //channel_info[chnum].ReadWf(waveform.substr(waveform.find(" ")+1));

    // ////////##########/////////Starting analysis
    //         short_channel_info[chnum].Initialize();
                channel_info[chnum].OnlyWf(Event_Length);
                Int_t zero_level = channel_info[chnum].Get_Zero_Level(GATE_BEG);
                cout << zero_level << endl << endl << "wefwewefwsdvegbvsdvwdv" << endl; 
                channel_info[chnum].ShiftWf(Event_Length,zero_level);
                zero_level = 0;
            //  cout << " wejfnwie" << endl;
            // for (int i = 0; i < Event_Length; i++) cout << channel_info[0].wf[i] << " :: ";
            // cout << endl << endl << endl;               
            //channel_info[chnum].SplineWf();
            short_channel_info[chnum].zl_rms = channel_info[chnum].Get_Zero_Level_RMS(GATE_BEG);
            short_channel_info[chnum].integral_in_gate = channel_info[chnum].Get_Integral_in_Gate(zero_level, GATE_BEG, GATE_END);
        

            short_channel_info[chnum].peak_pos = 
            channel_info[chnum].Get_peak_position_gauss(zero_level, short_channel_info[chnum].amp, GATE_END);
            short_channel_info[chnum].amp = channel_info[chnum].Get_Amplitude(zero_level, GATE_BEG, GATE_END);//inv_amp;
            //channel_info[chnum].SplineWf();channel_info[chnum].SplineWf();channel_info[chnum].DiffWf(); channel_info[chnum].FindDiffWfPars(min_diff,min_time,max_diff,max_time); if (min_diff > -15 && max_diff < 15) short_channel_info[chnum].Initialize();


            combined_tree->Fill();
            event_number++;

            for (Int_t i = 0; i < total_channels; i++) channel_info[i].Initialize();
    }
combined_tree->Write();
combined_root->Close();
}

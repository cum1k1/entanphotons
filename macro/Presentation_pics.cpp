#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>

#include <TF1.h>
#include <TCanvas.h>
#include <RooDouble.h>

#include "like_ivashkin_wants_it.h"

#include <TROOT.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLeaf.h>

#include <TChain.h>
#include <TString.h>
#include <TDirectory.h>
#include <TSystemFile.h>
#include <TSystemDirectory.h>
#include <TGraphErrors.h>
#include <TObjString.h>
#include <TPaveText.h>
#include "ResetHistogramPeakToZero.h"
#include "ChannelEntry.h"
#include "Check_entanglement_cuts.h"

using namespace HistPeakToZero;

void Presentation_pics()
{
    TString source_path = "/home/doc/entanglement/with_spline/decoherent/";
    TH1F *shifted_nrg = new TH1F("shifted_nrg","",1,0,1);
    hist_like_ivashkin_wants_it(shifted_nrg,"Time [ns]", "Normalized Number of Events", 1);
    const Int_t total_channels = 35;
    TChain *PMT_tree = new TChain;
	PMT_tree->AddFile( (source_path + "calibrated_time.root" + "/adc64_data").Data() );
    TFile *entangled_file = TFile::Open(source_path + "entangled_mini_tree.root","readonly");
	auto *entangled_tree = entangled_file->Get<TTree>("Signals");    
    TFile *merged_file = TFile::Open("~/Downloads/merged.root","readonly");
	auto *merged_tree = merged_file->Get<TTree>("Signals");

    TFile *decoherent_file = TFile::Open(source_path + "decoh_mini_tree.root","readonly");
	auto *decoherent_tree = entangled_file->Get<TTree>("Signals");

    TFile *data_file = TFile::Open(source_path+"07a8de9a_20211105_014338.root","readonly");
    auto * data_tree = data_file->Get<TTree>("adc64_data");

    Double_t factor = 1;
    Int_t samples[2050] = {0};            Int_t initial_point = 40;
    //TF1 * gaus_fit = new TF1("gaus_fit","gaus",-80,80); 
    Float_t low_cut, high_cut;

    ChannelEntry *channel_info = new ChannelEntry[total_channels];
    for(Int_t ch = 0; ch < total_channels; ch++)
	    (channel_info+ch)->SetBranch(data_tree, ch);

    for (Int_t i = 0; i < 350; i++) samples[i] = i*16;

    Int_t draw_events = 200000000;

///////////////###########Drawing GAGG
    TCanvas *gagg_canv = new TCanvas("gagg_canv", "gagg_canv");
    TH1F *gagg_hist = new TH1F("gagg_hist","", 200,-100,100);
    hist_like_ivashkin_wants_it(gagg_hist,"Time [ns]", "Normalized Number of Events", 1);
    PMT_tree->Draw("channel_34.peak_pos - channel_32.peak_pos >> gagg_hist","channel_32.integral_in_gate > 100 && channel_32.integral_in_gate < 350"
    "&& channel_34.integral_in_gate > 4 && channel_34.integral_in_gate < 50"
    "&& diff_channel_34.min_diff < -7 && diff_channel_34.max_diff > 7","",draw_events);
    gagg_hist->Scale(factor/gagg_hist->GetBinContent(gagg_hist->GetMaximumBin()));
    TH1F *gagg_hist_be = new TH1F("gagg_hist_be","gagg_hist_be", 200,-50,150);
    hist_like_ivashkin_wants_it(gagg_hist_be,"Time [ns]", "Normalized Number of Events", 2);
    PMT_tree->Draw("channel_34.peak_pos - channel_32.peak_pos >> gagg_hist_be","channel_32.integral_in_gate > 100 && channel_32.integral_in_gate < 350"
    "&& channel_34.integral_in_gate > 50 && channel_34.integral_in_gate < 140"

    "&& diff_channel_34.min_diff < -7 && diff_channel_34.max_diff > 7","",draw_events);
    gagg_hist_be->Scale(factor/gagg_hist_be->GetBinContent(gagg_hist_be->GetMaximumBin()));
    *shifted_nrg = ShiftedHist(gagg_hist_be,gagg_hist_be->GetMaximumBin()-find_zero_bin(gagg_hist_be),3,-40,40);
    TH1F *shifted_else = new TH1F("shifted_else","shifted_else",1,0,1);
    *shifted_else = ShiftedHist(gagg_hist,gagg_hist_be->GetMaximumBin()-find_zero_bin(gagg_hist_be),3,-40,40);

    hist_like_ivashkin_wants_it(shifted_nrg,"Time [ns]", "Normalized Number of Events", 1);
    hist_like_ivashkin_wants_it(shifted_else,"Time [ns]", "Normalized Number of Events", 2);
    TF1 *gaus_fit = new TF1("gaus_fit","gaus",-5.7,6);
    gaus_fit->SetLineColor(3);
    TF1 *gaus_fit_2 = new TF1("gaus_fit_2","gaus",-8,12.5);
    gaus_fit_2->SetLineColor(4);

    shifted_nrg->Fit("gaus_fit","R");
    shifted_else->Fit("gaus_fit_2","R");
    
    shifted_nrg->Draw("");
    shifted_else->Draw("same");
    gagg_canv->Update();
    gagg_canv->SaveAs(source_path+"cool_pics.pdf(",".pdf");

    delete gagg_hist;
    delete gagg_hist_be;
    delete gagg_canv;
    delete gaus_fit;
    ///###???
    /////Drawing main sc
    gagg_canv = new TCanvas("gagg_canv", "gagg_canv");
    gagg_hist = new TH1F("gagg_hist","", 350,0,700);
    hist_like_ivashkin_wants_it(gagg_hist,"Energy [keV]", "Normalized Number of Events", 1);
    PMT_tree->Draw("channel_32.integral_in_gate >> gagg_hist","abs(channel_32.peak_pos - channel_33.peak_pos) < 7"
    "&& channel_32.peak_pos > 10 && channel_33.peak_pos > 10 && channel_32.amp > 1000","",draw_events);
    gagg_hist->Scale(factor/gagg_hist->GetBinContent(gagg_hist->GetMaximumBin()));

    gagg_hist_be = new TH1F("gagg_hist_be","gagg_hist_be", 350,0,700);
    hist_like_ivashkin_wants_it(gagg_hist_be,"Energy [keV]", "Normalized Number of Events", 2);
    entangled_tree->Draw("EdepScat0 >> gagg_hist_be","EdepDet0 > 200 && EdepDet0 < 300"
    ,"",draw_events);
    gagg_hist_be->Scale(factor/gagg_hist_be->GetBinContent(gagg_hist_be->GetMaximumBin()));

    gagg_hist->Draw("hist");
    gagg_hist_be->Draw("histsame");
    gagg_canv->Update();
    gagg_canv->SaveAs(source_path+"cool_pics.pdf",".pdf");

    delete gagg_hist;
    delete gagg_hist_be;
    delete gagg_canv;
    ///###???   

    /////Drawing GAGG and noise
    gagg_canv = new TCanvas("gagg_canv", "gagg_canv");
    //gagg_hist = new TH1F("gagg_hist","gagg_hist", 350,0,700);
    decoherent_file->GetObject("integral_in_gate_diffuser",gagg_hist);
    hist_like_ivashkin_wants_it(gagg_hist,"Energy [keV]", "Normalized Number of Events", 1);
    gagg_hist->Scale(2*factor/(gagg_hist->GetBinContent(1)+gagg_hist->GetBinContent(2)+gagg_hist->GetBinContent(3)));
    gagg_hist->SetTitle("");
    entangled_file->GetObject("integral_in_gate_diffuser",gagg_hist_be);
    hist_like_ivashkin_wants_it(gagg_hist_be,"Energy [keV]", "Normalized Number of Events", 2);
    gagg_hist_be->Scale(2*factor/(gagg_hist_be->GetBinContent(1)+gagg_hist_be->GetBinContent(2)+gagg_hist_be->GetBinContent(3)));
    gagg_hist->GetXaxis()->SetRangeUser(0,140);
    gagg_hist_be->GetXaxis()->SetRangeUser(0,140);
    gagg_hist->Draw("hist");
    gagg_hist_be->Draw("histsame");
    gagg_canv->Update();
    gagg_canv->SaveAs(source_path+"cool_pics.pdf",".pdf");

    delete gagg_hist;
    delete gagg_hist_be;
    delete gagg_canv;
   ////////////////////////////
    //NaI spectrum
    TCanvas *NaI_canv = new TCanvas("NaI_canv", "NaI_canv");
    TH1F *NaI_hist = new TH1F("NaI_hist","", 200,0,400);
    TH1F *mcNaI_hist = new TH1F("mcNaI_hist","",200,0,400);
    hist_like_ivashkin_wants_it(NaI_hist,"Energy [keV]", "Normalized Number of Events", 1);
    hist_like_ivashkin_wants_it(mcNaI_hist,"Energy [keV]", "Normalized Number of Events", 2);
   
    entangled_tree->Draw("EdepDet0 >> NaI_hist","EdepScat0 > 200 && EdepScat0 < 300");
    NaI_hist->Scale(factor/NaI_hist->GetEntries());

    merged_tree->Draw("MiniTree.EdepDet0 >> mcNaI_hist");
    mcNaI_hist->Scale(factor/mcNaI_hist->GetEntries());
    Float_t scaling = NaI_hist->GetBinContent(mcNaI_hist->GetMaximumBin());
    // mcNaI_hist->Scale(1./scaling);
     NaI_hist->Scale(1./scaling);
    // mcNaI_hist->Draw("hist");

    NaI_hist->Draw("hist");
    NaI_canv->Update();
    NaI_canv->SaveAs(source_path+"cool_pics.pdf",".pdf");
//////////////////////////
    //NaI_vs_sc spectrum
    TCanvas *NaI_canv2d = new TCanvas("NaI_canv2d", "");
    TH2F *NaI_vs_sc_hist = new TH2F("NaI_vs_sc_hist","", 200,180,350,200,200,320);
    hist2d_like_ivashkin_wants_it(NaI_vs_sc_hist,"Energy in scatterer [keV]", "Energy in NaI [keV]", 1);
    entangled_tree->Draw("EdepDet0:EdepScat0 >> NaI_vs_sc_hist","","colz");
    NaI_canv2d->SetLogz();
    //NaI_vs_sc_hist->Draw("hist");
    NaI_canv2d->Update();
    NaI_canv2d->SaveAs(source_path+"cool_pics.pdf",".pdf");

    delete NaI_vs_sc_hist;
    delete NaI_canv2d;

//////////////////////////
//##32-33 difference time
    gagg_canv = new TCanvas("gagg_canv", "gagg_canv");
    gagg_hist = new TH1F("gagg_hist","", 21,-10,10);
    hist_like_ivashkin_wants_it(gagg_hist,"Time [ns]", "Normalized Number of Events", 1);
    PMT_tree->Draw("channel_33.peak_pos - channel_32.peak_pos >> gagg_hist","channel_2.integral_in_gate > 100","",draw_events);
    gagg_hist->Scale(factor/gagg_hist->GetBinContent(gagg_hist->GetMaximumBin()));
    *shifted_nrg = ShiftedHist(gagg_hist,gagg_hist->GetMaximumBin()-find_zero_bin(gagg_hist));
    shifted_nrg->GetXaxis()->SetRangeUser(-10,10);
    shifted_nrg->GetYaxis()->SetRangeUser(0,1.1);
    hist_like_ivashkin_wants_it(shifted_nrg,"Time [ns]", "Normalized Number of Events", 1);
    gaus_fit = new TF1("gaus_fit","gaus",-1.8,3.2);
    gaus_fit->SetLineColor(3);
    shifted_nrg->Fit("gaus_fit","R");
    shifted_nrg->Draw("");
    gagg_canv->Update();
    gagg_canv->SaveAs(source_path+"cool_pics.pdf",".pdf");
    /////////////
    //##32-2 difference time
    gagg_canv = new TCanvas("gagg_canv", "gagg_canv");
    gagg_hist = new TH1F("gagg_hist","", 901,-450,450);
    hist_like_ivashkin_wants_it(gagg_hist,"Time [ns]", "Normalized Number of Events", 1);
    PMT_tree->Draw("channel_2.peak_pos - channel_32.peak_pos >> gagg_hist","channel_32.integral_in_gate > 200 && channel_32.integral_in_gate < 300"
    "&& channel_2.integral_in_gate > 200 && channel_2.integral_in_gate < 300","",draw_events);
    gagg_hist->Scale(factor/gagg_hist->GetBinContent(gagg_hist->GetMaximumBin()));
    *shifted_nrg = ShiftedHist(gagg_hist,gagg_hist->GetMaximumBin()-find_zero_bin(gagg_hist));
    hist_like_ivashkin_wants_it(shifted_nrg,"Time [ns]", "Normalized Number of Events", 1);
    cout << find_zero_bin(gagg_hist) << " " << gagg_hist->FindBin(0);
    gaus_fit = new TF1("gaus_fit","gaus",-15,15);
    gaus_fit->SetLineColor(3);
    shifted_nrg->Fit("gaus_fit","R");

    shifted_nrg->Draw("");
    TPaveText *pt = new TPaveText(15.35,0.5,28.95,0.6);
    pt->AddText(Form("#sigma_{t} = %4.3f",gaus_fit->GetParameter(2)));    
    //pt->Draw("same");
    gagg_canv->Update();
    gagg_canv->SaveAs(source_path+"cool_pics.pdf",".pdf");
    /////////////

    ////////////Drawing waveforms
        Int_t global_counter = 0;
        Int_t tot_canv = 1;
        for (Int_t entryNum = 0; entryNum < data_tree->GetEntries()/10; entryNum++)
        {
            data_tree->GetEntry(entryNum);

            const Int_t wfs_34 = channel_info[34].wf_size; Int_t wf_34[2000] = {0};
            const Int_t wfs_32 = channel_info[32].wf_size; Int_t wf_32[2000] = {0};
            const Int_t wfs_5 = channel_info[5].wf_size; Int_t wf_5[2000] = {0};
            Int_t amp34 = channel_info[34].Get_Amplitude(channel_info[34].Get_Zero_Level(50), 50, 85);
            Int_t amp32 = channel_info[32].Get_Amplitude(channel_info[32].Get_Zero_Level(50), 50);
            Int_t amp5 = channel_info[5].Get_Amplitude(channel_info[5].Get_Zero_Level(50), 50,150);                        


            //for (Int_t i = initial_point; i <= wfs_5; i++)    (Int_t)(wf_5 [i-initial_point] = channel_info[ 5].Get_Zero_Level(50)-channel_info[ 5].wf[i]);
                    //for (Int_t i = 0; i < 150; i++) wf_5[i] = (Int_t)channel_info[5].Get_Zero_Level(50)-(Int_t)channel_info[5].wf[i];

            Short_t flag = 0;
            if (global_counter < tot_canv)
            {
                if (wfs_34 > 10 
                && wfs_32 > 10
                && wfs_5 > 10
                && amp34 > 2000
                && amp32 > 5000
                && amp5 > 7000
                && amp34 < 10000
                && amp32 < 30000
                && amp5 < 30000                )
                {
                    
                    if (global_counter < tot_canv)
                    {
                        TCanvas *canv = new TCanvas("canv","canv");

                        for (Int_t i = initial_point; i <= wfs_34; i++) wf_34[i-initial_point] = channel_info[34].Get_Zero_Level(50)-channel_info[34].wf[i];
                        for (Int_t i = initial_point; i <= 71; i++) wf_32[i-initial_point] = channel_info[32].Get_Zero_Level(50)-channel_info[32].wf[i];for (Int_t i = 71-initial_point; i < 2000; i++) wf_32[i]=0;
                        for (Int_t i = initial_point; i <= wfs_5; i++) wf_5[i-initial_point] = channel_info[5].Get_Zero_Level(50)-channel_info[5].wf[i];
                        TGraph *gr_32 = new TGraph(150-40,samples,wf_32);
                        graph_like_ivashkin_wants_it(gr_32,"Time [ns]","ADC channels", "",1,1);
                            gr_32->GetXaxis()->SetRangeUser(0,1200);

                        gr_32->Draw("AL");        

                        TGraph *gr_5 = new TGraph(150-40,samples,wf_5);
                        graph_like_ivashkin_wants_it(gr_5,"Time [ns]","ADC channels","",8,8);
                        gr_5->GetXaxis()->SetRangeUser(0,1200);
                        gr_5->Draw("L");     

                        TGraph *gr_34 = new TGraph(150-40,samples,wf_34);
                        graph_like_ivashkin_wants_it(gr_34,"Time [ns]","ADC channels", "",2,2);
                        gr_34->GetXaxis()->SetRangeUser(0,1200);
                        gr_34->Draw("L");    
                        canv->SaveAs(source_path+"cool_pics.pdf","pdf");
                    delete canv;
                        global_counter++;
                    }    

                }
            }
            else
            {entryNum = data_tree->GetEntries()-1;}
            }

    //////////////
    /////Drawing GAGG and noise
    gagg_canv = new TCanvas("gagg_canv", "gagg_canv");
    gagg_hist = new TH1F("gagg_hist","gagg_hist", 290,10,300);
    PMT_tree->Draw("channel_34.integral_in_gate >> gagg_hist","channel_34.peak_pos - channel_32.peak_pos < 80 && channel_34.peak_pos - channel_32.peak_pos > 40","",2*draw_events);
    hist_like_ivashkin_wants_it(gagg_hist,"Energy [keV]", "Normalized Number of Events", 1);
    gagg_hist->Scale(1/(gagg_hist->GetBinContent(gagg_hist->GetMaximumBin())));
    gagg_hist->SetTitle("");
    gagg_hist->GetXaxis()->SetRangeUser(10,300);
    gagg_hist->Draw("hist");
    gagg_canv->Update();
    gagg_canv->SaveAs(source_path+"cool_pics.pdf)",".pdf");

    delete gagg_hist;
    delete gagg_hist_be;
    delete gagg_canv;
   ////////////////////////////
}
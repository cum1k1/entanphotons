// #include <iostream>
// #include <fstream>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCutG.h>

// #include <TObjArray.h>
#include <TString.h>
#include <TGraphErrors.h>
#include "ChannelEntry.h"
#include "CHSH_calculator.h"
#include "Check_entanglement_cuts.h"
#include "like_ivashkin_wants_it.h"
using namespace std;
using namespace CHSH;
using namespace CUTS;

#define UseDecoherentPhotons 1
#define Analyze_Monte_Carlo 0
#define UseNoCutG 0


void MiniDST_analysis()
{
	TString source_path = "/home/doc/entanglement/with_spline/decoherent/";

#if UseDecoherentPhotons
    TString middle_path  = source_path + "decoh_mini_tree";
#else
    TString middle_path = source_path + "entangled_mini_tree";
#endif
    mini_tree_nrg *short_tree = new mini_tree_nrg;
    TFile *f = TFile::Open(middle_path+".root");
    TTree *MiniDST_tree = (TTree*)f->Get("Signals");
    
    TString result_path = middle_path + "_analyzed";
    TFile *result_root = new TFile(result_path+".root","RECREATE");
    short_tree->SetBranch(MiniDST_tree);

    Float_t total_events_for_angle_diff[16] = {0};
    Float_t total_events_for_angle_diff_err[16] = {0};    
    Int_t NumEvents[16][16] = {0};     
    Float_t angle_arr[16]  = {0.};
    Float_t angle_arr_err[16] = {0.};    

    for (Int_t i = 0; i <16; i++) angle_arr[i] = (float)i*22.5;
    /////////////////////////////
    //////////////////////////////

        Float_t low_det0_cut = 180;
        Float_t high_det0_cut = 300;
        Float_t low_det1_cut = 180;
        Float_t high_det1_cut = 300;

        Float_t low_scat0_cut = 150;
        Float_t high_scat0_cut = 300; 
        Float_t low_scat1_cut = 180;
        Float_t high_scat1_cut  = 300;

        Float_t high_Intermediate_cut = 60;
        Float_t low_Intermediate_cut = 2;
/////////////////////////////////
////////////////////////////////

//////////////Setting TCUTG boarders
    TCanvas *canv_0 = new TCanvas("canv_0,canv_0");
    canv_0->cd();

    #if UseDecoherentPhotons
    //MiniDST_tree->Draw("MiniTree.EdepDet0:MiniTree.EdepIntermediate >> (100,0,300,100,0,400)"/*,Form("MiniTree.EdepIntermediate < %f && MiniTree.EdepIntermediate > %f", high_Intermediate_cut, low_Intermediate_cut)*/,"", "colz");
        MiniDST_tree->Draw("MiniTree.EdepDet0:MiniTree.EdepScat0 >> (75,0,500,75,0,500)","", "colz");

    #else
    MiniDST_tree->Draw("MiniTree.EdepDet0:MiniTree.EdepScat0 >> (300,0,500,300,0,500)","", "colz");
    #endif
    
   TCutG *cutg = new TCutG("CUTG",13);
   cutg->SetVarX("EdepScat0 ");
   cutg->SetVarY("EdepDet0");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->SetPoint(0,156.148,277.75);
   cutg->SetPoint(1,181.596,241.838);
   cutg->SetPoint(2,211.116,208.829);
   cutg->SetPoint(3,241.653,193.594);
   cutg->SetPoint(4,275.448,188.516);
   cutg->SetPoint(5,301.303,219.711);
   cutg->SetPoint(6,287.663,255.622);
   cutg->SetPoint(7,248.168,279.201);
   cutg->SetPoint(8,214.78,301.328);
   cutg->SetPoint(9,182.003,327.445);
   cutg->SetPoint(10,150.041,332.886);
   cutg->SetPoint(11,142.915,303.867);
   cutg->SetPoint(12,156.148,277.75);
   cutg->Draw("same");
    canv_0->SaveAs(result_path+".pdf(",".pdf");
    canv_0->Write();

//////////////////////////////////////////
        
        TH1F *hist_for_fit = new TH1F("hist_for_fit","hist_for_fit", 200, 50, 450);
        TCanvas *new_canv = new TCanvas("new", "new");
        MiniDST_tree->Draw("MiniTree.EdepDet0 >> hist_for_fit"

#if UseDecoherentPhotons

        ,Form("MiniTree.EdepIntermediate < %f && MiniTree.EdepIntermediate > %f", high_Intermediate_cut, low_Intermediate_cut)
#endif        
        ); 
        double_gauss_fit(hist_for_fit, low_det0_cut, high_det0_cut);
        hist_for_fit->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();

        MiniDST_tree->Draw("MiniTree.EdepDet1 >> hist_for_fit");
        double_gauss_fit(hist_for_fit, low_det1_cut, high_det1_cut);
        hist_for_fit->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();

        MiniDST_tree->Draw("MiniTree.EdepScat1 >> hist_for_fit", 
        Form("MiniTree.EdepDet1 < %f && MiniTree.EdepDet1 > %f", high_det1_cut, low_det1_cut));
        double_gauss_fit(hist_for_fit, low_scat1_cut, high_scat1_cut);
        hist_for_fit->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();


        MiniDST_tree->Draw("MiniTree.EdepScat0 >> hist_for_fit", 
        Form("MiniTree.EdepDet0 < %f && MiniTree.EdepDet0 > %f" 
#if UseDecoherentPhotons
        "&& MiniTree.EdepIntermediate < %f && MiniTree.EdepIntermediate > %f"
#endif
        , high_det0_cut, low_det0_cut 
#if UseDecoherentPhotons
        ,high_Intermediate_cut, low_Intermediate_cut
#endif
        ));
        double_gauss_fit(hist_for_fit, low_scat0_cut, high_scat0_cut);
        hist_for_fit->Write();
        new_canv->SaveAs(result_path+".pdf(",".pdf");
        new_canv->Clear();

   


/////////////////////Selecting coicidences for all counters
    for (Int_t NumEvent = 0; NumEvent < MiniDST_tree->GetEntries(); NumEvent++)
    {
        MiniDST_tree->GetEntry(NumEvent);     
        Short_t num0 = short_tree->DetNum0;
        Short_t num1 = short_tree->DetNum1;
        cout << num0 << " "<<num1<<endl;

        #if Analyze_Monte_Carlo
        num0-- ; num1--;
        #endif
        if(short_tree->EdepDet1 > low_det1_cut
        && short_tree->EdepDet1 < high_det1_cut
        && short_tree->EdepScat1 > low_scat1_cut
        && short_tree->EdepScat1 < high_scat1_cut 
#if UseNoCutG
        && short_tree->EdepDet0 > low_det0_cut
        && short_tree->EdepDet0 < high_det0_cut        
        && short_tree->EdepScat0 > low_scat0_cut
        && short_tree->EdepScat0 < high_scat0_cut
#else
        //&& short_tree->EdepIntermediate > 10
        //&& short_tree->EdepDet0 < 340
        //&& short_tree->EdepDet0 > 270 
        //&& short_tree->EdepIntermediate > 25
        //&& short_tree->EdepIntermediate < 125

        && cutg->IsInside(short_tree->EdepScat0,short_tree->EdepDet0)
        
#endif 

#if UseDecoherentPhotons
        && short_tree->EdepIntermediate < high_Intermediate_cut
        && short_tree->EdepIntermediate > low_Intermediate_cut   
#endif
        && num0 > -1 && num1 > -1)
            NumEvents[num0][num1-16] +=1;
    }
    //#############end of actual program further only CHSH and ratio calculation and drawing
    /////////////////////
    /////////////////////Calculate coincidences for ratio

        for (Int_t channel_number = 0; channel_number < 16; channel_number++)
        {
            for(Int_t channel_number_2 = 16; channel_number_2 < 32; channel_number_2++)
            {
                Int_t angle = channel_number_2-channel_number-16;
                if (angle < 0) angle += 16;    
                total_events_for_angle_diff[angle]+=NumEvents[channel_number][channel_number_2-16];
                // Int_t i = true_number(channel_number);
                // Int_t j = true_number(channel_number_2);
                // Int_t i1 = true_number(channel_number-2);
                // Int_t j1 = true_number(channel_number_2-2);
                // total_events_for_angle_diff[angle]+=NumEvents[i][j]+NumEvents[j][i]+NumEvents[i][j1]+NumEvents[j][i1];
            
            }
        }

/////////////////////
/////////////////////CHSH inequality and Correlation coefficients
        TF1 *my_fit = new TF1 ("my_fit","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",-190,190);
        TF1 *e_fit = new TF1("e_fit","[0]*(cos(2*3.141593/180*x))",-190,190);

            const Int_t total_angles = 8;
            const Int_t all_angles = total_angles*2;
            Float_t phi_angle[all_angles] = {0.};
            for (Int_t i = 0; i < total_angles; i++) phi_angle[i] = (float)(i+1)*22.5;
            for (Int_t i = 0; i < total_angles; i++) phi_angle[total_angles+i] = -(float)(i+1)*22.5;

            Float_t CHSH[total_angles] = {0}; Float_t CHSH_all[all_angles] = {0};
            Float_t CHSH_error[total_angles] = {0}; Float_t CHSH_error_all[all_angles] = {0};
            Float_t global_CHSH[total_angles] = {0}; Float_t global_CHSH_all[all_angles] = {0};
            Float_t global_CHSH_error[total_angles] = {0}; Float_t global_CHSH_error_all[all_angles] = {0};
            Float_t global_Correlation_E[total_angles] = {0}; Float_t global_Corr_error[total_angles] = {0};
            Float_t global_Correlation_E_all[all_angles] = {0}; Float_t global_Corr_error_all[all_angles] = {0};

            for (Int_t a_angle = 0; a_angle < 16; a_angle++)
            {
                Float_t CHSH_error_local[all_angles] = {0};  Float_t CHSH_error_local_sum[total_angles] = {0};
                Float_t CHSH_local[all_angles] = {0}; Float_t CHSH_local_sum[total_angles] = {0};
                for (Int_t angle_between_counters = 1; angle_between_counters <= total_angles; angle_between_counters++)
                {
                    CHSH_local[angle_between_counters-1] = 
                    calculate_local_CHSH(NumEvents, a_angle, angle_between_counters, "clockwise");
                    CHSH_error_local[angle_between_counters-1] = 
                    calculate_local_CHSH_error(NumEvents, a_angle, angle_between_counters, "clockwise");
                    CHSH_local[total_angles + angle_between_counters-1] = 
                    calculate_local_CHSH(NumEvents, a_angle, angle_between_counters, "counterclockwise");
                    CHSH_error_local[total_angles + angle_between_counters-1] = 
                    calculate_local_CHSH_error(NumEvents, a_angle, angle_between_counters, "counterclockwise");

                    CHSH_local_sum[angle_between_counters-1] = (CHSH_local[angle_between_counters-1] + CHSH_local[angle_between_counters+total_angles-1])/2;
                    CHSH_error_local_sum[angle_between_counters-1] = CHSH_error_local[angle_between_counters-1] + CHSH_error_local[angle_between_counters+total_angles-1];
                    
                    CHSH[angle_between_counters-1] += CHSH_local_sum[angle_between_counters-1];
                    CHSH_error[angle_between_counters-1] += CHSH_error_local_sum[angle_between_counters-1];       
                    CHSH_error_local_sum[angle_between_counters-1] = sqrt(CHSH_error_local_sum[angle_between_counters-1])/2;  

                    CHSH_all[angle_between_counters-1] += CHSH_local[angle_between_counters-1];
                    CHSH_all[angle_between_counters + total_angles-1] += CHSH_local[angle_between_counters + total_angles-1];

                    CHSH_error_all[angle_between_counters-1] += CHSH_error_local[angle_between_counters-1];  
                    CHSH_error_all[angle_between_counters + total_angles-1] += CHSH_error_local[angle_between_counters + total_angles-1];  

                    CHSH_error_local[angle_between_counters-1] = sqrt(CHSH_error_local[angle_between_counters-1]); 
                    CHSH_error_local[angle_between_counters + total_angles-1] = sqrt(CHSH_error_local[angle_between_counters + total_angles-1]); 

                    global_Correlation_E_all[angle_between_counters-1] += E_coeff(NumEvents, a_angle, a_angle + angle_between_counters);
                    global_Corr_error_all[angle_between_counters-1] += sqr_E_error(NumEvents, a_angle, a_angle + angle_between_counters);
                    global_Correlation_E_all[angle_between_counters+total_angles-1] += E_coeff(NumEvents, a_angle, a_angle + angle_between_counters, "counterclockwise");
                    global_Corr_error_all[angle_between_counters+total_angles-1] += sqr_E_error(NumEvents, a_angle, a_angle + angle_between_counters, "counterclockwise");
                }

                    TCanvas *canvas_CHSH = new TCanvas ( Form("canvas_CHSH_a=%i",a_angle), Form("canvas_CHSH_a=%i",a_angle));
                    TGraphErrors * gr_CHSH = new TGraphErrors(total_angles,phi_angle,CHSH_local_sum, angle_arr_err, CHSH_error_local_sum);
                    //TF1 *my_fit_2 = new TF1 ("my_fit_2","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",10,190);
                    canvas_CHSH->Divide(2);

                    gr_CHSH->Fit("my_fit","R");
                    Float_t a_par_CHSH = abs(my_fit->GetParameter(0));
                    Float_t a_error_CHSH = abs(my_fit->GetParError(0));
                    graph_like_ivashkin_wants_it(gr_CHSH,"angle [degrees]","S", 
                    Form("canvas_CHSH_a=%i <> E(90) = %4.3f+-%4.3f",a_angle,
                    (E_coeff(NumEvents,a_angle,true_number(a_angle+4))+E_coeff(NumEvents,a_angle,true_number(a_angle-4)))/2,
                    sqrt(sqr_E_error(NumEvents,a_angle,true_number(a_angle+4))+sqr_E_error(NumEvents,a_angle,true_number(a_angle-4)))/2),1);
                    canvas_CHSH->cd(1);
                    gr_CHSH->Draw("AP");

                    //TF1 *my_fit = new TF1 ("my_fit","[0]*(cos(6*3.141593/180*x)-3*cos(2*3.141593/180*x))",-190,190);
                    TGraphErrors * gr_CHSH_all = new TGraphErrors(all_angles-1,phi_angle,CHSH_local, angle_arr_err, CHSH_error_local);
                    gr_CHSH_all->Fit("my_fit","R");
                    a_par_CHSH = abs(my_fit->GetParameter(0));
                    a_error_CHSH = abs(my_fit->GetParError(0));
                    graph_like_ivashkin_wants_it(gr_CHSH_all,"angle [degrees]","S", Form("canvas_CHSH_a=%i <> E(90) = %4.3f+-%4.3f",a_angle,E_coeff(NumEvents,a_angle,true_number(a_angle+4)),sqrt(sqr_E_error(NumEvents,a_angle,true_number(a_angle+4)))),1);
                    canvas_CHSH->cd(2);
                    gr_CHSH_all->Draw("AP");
                    gr_CHSH_all->Write();
                    canvas_CHSH->SaveAs(result_path+".pdf",".pdf");
                    delete canvas_CHSH;
            }

            cout <<endl<<endl<<endl<<endl;
            for (Int_t i = 0; i < all_angles; i++) 
            {

                if ( i < total_angles)
                {
                    CHSH[i] = CHSH[i]/16;
                    CHSH_error[i] = sqrt(CHSH_error[i])/32;
                    global_Correlation_E[i] = (global_Correlation_E_all[i] + global_Correlation_E_all[i+total_angles])/32;
                    global_Corr_error[i] = sqrt(global_Corr_error_all[i]+ global_Corr_error_all[i+total_angles])/32;
                }
                    CHSH_all[i] = CHSH_all[i]/16;
                    CHSH_error_all[i] = sqrt(CHSH_error_all[i])/16;
                    global_Corr_error_all[i] = sqrt(global_Corr_error_all[i])/16;
                    global_Correlation_E_all[i] = global_Correlation_E_all[i]/16;

            }
            for (Int_t i = 0; i < 8; i++) cout<< CHSH[i] <<" ";
            cout <<endl<<endl<<endl<<endl; 
            for (Int_t i = 0; i < 8; i++) cout<< CHSH_error[i] <<" ";
            cout <<endl<<endl<<endl<<endl;                 

            cout <<endl<<endl<<endl<<endl;

            Float_t total_Correlation_E[total_angles] = {0}; Float_t total_Corr_error[total_angles] = {0};
            Float_t total_Correlation_E_all[all_angles] = {0}; Float_t total_Corr_error_all[all_angles] = {0};

            for (Int_t angle_between_counters = 1; angle_between_counters <= total_angles; angle_between_counters++)
            {
                global_CHSH_all[angle_between_counters-1] = global_calculate_CHSH(angle_between_counters, NumEvents,"clockwise");
                global_CHSH_error_all[angle_between_counters-1] = sqrt(global_calculate_CHSH_error(angle_between_counters, NumEvents,"clockwise"));
                global_CHSH_all[total_angles+angle_between_counters-1] = global_calculate_CHSH(angle_between_counters, NumEvents,"counterclockwise");
                global_CHSH_error_all[total_angles+angle_between_counters-1] = sqrt(global_calculate_CHSH_error(angle_between_counters, NumEvents,"counterclockwise"));

                total_Correlation_E_all[angle_between_counters-1] = global_E_coeff(angle_between_counters, NumEvents,"clockwise");
                total_Corr_error_all[angle_between_counters-1] = sqrt(global_sqr_E_error(angle_between_counters, NumEvents,"clockwise"));
                total_Correlation_E_all[total_angles +angle_between_counters-1] = global_E_coeff(angle_between_counters, NumEvents,"counterclockwise");
                total_Corr_error_all[total_angles + angle_between_counters-1] = sqrt(global_sqr_E_error(angle_between_counters, NumEvents,"counterclockwise"));

                global_CHSH[angle_between_counters-1] = global_calculate_CHSH(angle_between_counters, NumEvents,"both");
                global_CHSH_error[angle_between_counters-1] = sqrt(global_calculate_CHSH_error(angle_between_counters, NumEvents,"both"));

                total_Correlation_E[angle_between_counters-1] = global_E_coeff(angle_between_counters, NumEvents,"both");
                total_Corr_error[angle_between_counters-1] = sqrt(global_sqr_E_error(angle_between_counters, NumEvents,"both"));

                cout<< global_CHSH_error[angle_between_counters-1] <<" ";
            }
            cout <<endl<<endl<<endl<<endl;
///////////////////////////////////////
///////////////////////////////////////
    Float_t average_E = 0; Float_t sigma_E_average = 0;

    for (Int_t a = 0; a < 16; a++)
    {
        Int_t b = true_number(a+4); Int_t b_1 = true_number(a-4);
        average_E += (E_coeff(NumEvents,a,b) + E_coeff(NumEvents,a,b_1))/32.;
        sigma_E_average += (sqr_E_error(NumEvents,a,b)+sqr_E_error(NumEvents,a,b_1));
    }
    sigma_E_average = sqrt(sigma_E_average)/32;
    TCanvas *canvas_CHSH = new TCanvas ( "canvas_CHSH", "canvas_CHSH");
    canvas_CHSH->Divide(2);

    TGraphErrors * gr_CHSH = new TGraphErrors(total_angles,phi_angle,CHSH, angle_arr_err, CHSH_error);
	gr_CHSH->Fit("my_fit","R");
    Float_t a_par_CHSH = abs(my_fit->GetParameter(0));
    Float_t a_error_CHSH = abs(my_fit->GetParError(0));
    canvas_CHSH->cd(1);
    graph_like_ivashkin_wants_it(gr_CHSH,"angle [degrees]","S", Form("Average E <> E(90) = %4.3f+-%4.3f",average_E,sigma_E_average),1);
    gr_CHSH->Draw("AP");
    TGraphErrors * gr_CHSH_all = new TGraphErrors(all_angles-1,phi_angle,CHSH_all, angle_arr_err, CHSH_error_all);
	gr_CHSH_all->Fit("my_fit","R");
    a_par_CHSH = abs(my_fit->GetParameter(0));
    a_error_CHSH = abs(my_fit->GetParError(0));
    canvas_CHSH->cd(2);
    graph_like_ivashkin_wants_it(gr_CHSH_all,"angle [degrees]","S", Form("Average E <> E(90) = %4.3f+-%4.3f",average_E,sigma_E_average),1);
    gr_CHSH_all->Draw("AP");
    gr_CHSH_all->Write();    
    canvas_CHSH->SaveAs(result_path+".pdf",".pdf");

    TCanvas *global_canvas_CHSH = new TCanvas ( "global_canvas_CHSH", "global_canvas_CHSH");
    global_canvas_CHSH->Divide(2);

    TGraphErrors * global_gr_CHSH = new TGraphErrors(total_angles,phi_angle,global_CHSH, angle_arr_err, global_CHSH_error);
	global_gr_CHSH->Fit("my_fit","R");
    a_par_CHSH = abs(my_fit->GetParameter(0));
    a_error_CHSH = abs(my_fit->GetParError(0));
    global_canvas_CHSH->cd(1);
    graph_like_ivashkin_wants_it(global_gr_CHSH,"angle [degrees]","S", Form("Sum N <> E(90) = %4.3f+-%4.3f",global_E_coeff(4,NumEvents,"both"), sqrt(global_sqr_E_error(4,NumEvents, "both"))),1);
    global_gr_CHSH->Draw("AP");
    
    TGraphErrors * global_gr_CHSH_all = new TGraphErrors(all_angles-1,phi_angle,global_CHSH_all, angle_arr_err, global_CHSH_error_all);
	global_gr_CHSH_all->Fit("my_fit","R");
    a_par_CHSH = abs(my_fit->GetParameter(0));
    a_error_CHSH = abs(my_fit->GetParError(0));
    global_canvas_CHSH->cd(2);
    graph_like_ivashkin_wants_it(global_gr_CHSH_all,"angle [degrees]","S", Form("Sum N <> E(90) = %4.3f+-%4.3f",global_E_coeff(4,NumEvents,"both"), sqrt(global_sqr_E_error(4,NumEvents, "both"))),1);
    
    global_gr_CHSH_all->Draw("AP");
    global_gr_CHSH_all->Write();
    global_canvas_CHSH->SaveAs(result_path+".pdf",".pdf");
  
   //////////////////////// Draw E_coeffs

    TCanvas *E_canv = new TCanvas ( "global_E", "global_E");
    E_canv->Divide(2);
    TGraphErrors * E_gr = new TGraphErrors(total_angles,phi_angle,global_Correlation_E, angle_arr_err, global_Corr_error);
    E_gr->Fit("e_fit","R");
    E_canv->cd(1);
    graph_like_ivashkin_wants_it(E_gr,"angle [degrees]","E", "Average E for mirroring angles", 1);
    E_gr->Draw("AP");

    TGraphErrors * E_gr_all = new TGraphErrors(all_angles-1,phi_angle,global_Correlation_E_all, angle_arr_err, global_Corr_error_all);
    E_gr_all->Fit("e_fit","R");
    E_canv->cd(2);
    graph_like_ivashkin_wants_it(E_gr_all,"angle [degrees]","E", "Average E for all angles", 1);
    E_gr_all->Draw("AP");
    E_gr_all->Write();    
    E_canv->SaveAs(result_path+".pdf",".pdf");

    TCanvas *E_total = new TCanvas ( "total_E", "total_E");
    E_total->Divide(2);
    TGraphErrors * E_gr_total = new TGraphErrors(total_angles,phi_angle, total_Correlation_E, angle_arr_err, total_Corr_error);
    E_gr_total->Fit("e_fit","R");
    E_total->cd(1);
    graph_like_ivashkin_wants_it(E_gr_total,"angle [degrees]","E", "Summ E for mirroring angles", 1);
    E_gr_total->Draw("AP");

    TGraphErrors * E_gr_all_total = new TGraphErrors(all_angles-1,phi_angle, total_Correlation_E_all, angle_arr_err, total_Corr_error_all);
    E_gr_all_total->Fit("e_fit","R");
    E_total->cd(2);
    graph_like_ivashkin_wants_it(E_gr_all_total,"angle [degrees]","E", "Summ E for all angles", 1);
    E_gr_all_total->Draw("AP");
    E_gr_all_total->Write();
    E_total->SaveAs(result_path+".pdf",".pdf");


   ///////////////////////Draw Ratio

    for (Int_t i = 0; i < 16; i++) total_events_for_angle_diff_err[i] = pow(total_events_for_angle_diff[i],0.5);

    TGraphErrors * gr = new TGraphErrors(16,angle_arr,total_events_for_angle_diff, angle_arr_err, total_events_for_angle_diff_err);
	TF1 *sin_fit = new TF1 ("sin_fit","[0]+[1]*cos(3.14159265359/180*2*x)",-10,360);
	gr->Fit("sin_fit","R");
    Float_t a_par = abs(sin_fit->GetParameter(0));
    Float_t b_par = abs(sin_fit->GetParameter(1));
    Float_t a_error = abs(sin_fit->GetParError(0));
    Float_t b_error = abs(sin_fit->GetParError(1));
    Float_t diff = (a_par+b_par)/(a_par-b_par);
    Float_t diff_err = 2*diff*sqrt(pow((a_error/a_par),2)+pow((b_error/b_par),2))/(a_par/b_par-b_par/a_par);
    //Float_t diff_err = 2*a_par/(pow(a_par-b_par,2))*sqrt(pow((a_error*b_par/a_par),2)+pow((b_error),2));

    TCanvas *canvas = new TCanvas ( "canvas", "canvas");
    canvas->cd();
    graph_like_ivashkin_wants_it(gr,"angle [degrees]","Counts", Form("%5.1f-%5.1f*cos(2x) ratio = %4.3f +- %4.3f",a_par,b_par,diff,diff_err),1);
    gr->Draw("AP");
    gr->Write();
    canvas->SaveAs(result_path+".pdf)",".pdf");    
    result_root->Close();    
    }




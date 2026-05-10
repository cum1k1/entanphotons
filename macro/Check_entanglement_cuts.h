#ifndef CHECK_CUTS
#define CHECK_CUTS
#include "ChannelEntry.h"

namespace CUTS
{

    Bool_t Decoherent_or_Entangled(short_energy_ChannelEntry *short_channel_info, Bool_t entangled = kTRUE, diff_short_energy_ChannelEntry *diff_short = nullptr)
    {
        Bool_t flag = kFALSE;
        if (entangled == kFALSE && short_channel_info[34].amp < 50000
            && short_channel_info[34].amp > 100
            && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 110
            && (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 20
            && diff_short->min_diff < -7 && diff_short->max_diff > 7
            )
            flag = kTRUE;

        if (entangled == kTRUE &&
            short_channel_info[34].amp < 300 &&
            ((short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) < 0 ||
            (short_channel_info[34].peak_pos - short_channel_info[32].peak_pos) > 150)
            )
            flag = kTRUE;
        return flag;
    }

    Bool_t Non_zero_time(short_energy_ChannelEntry *short_channel_info, Int_t channel_number)
    {
        Bool_t flag = kFALSE;
        if (
        short_channel_info[channel_number].peak_pos > 11
        &&  short_channel_info[33].peak_pos > 11
        &&  short_channel_info[32].peak_pos > 11) flag = kTRUE;
        return flag;
    }

    Bool_t Apply_time_in_peak_cuts(short_energy_ChannelEntry *short_channel_info, Int_t channel_number, Int_t sc_number = 32, TString data = "new")
    {
        Int_t main_sc_boarder = 0;
        Int_t left_ch_sc_boarder = 0;
        Int_t right_ch_sc_boarder = 0;
        if (data == "new")
        {
            main_sc_boarder = 9;
            left_ch_sc_boarder = 150;
            right_ch_sc_boarder = 300;
        }
        if (data == "old")
        {
            main_sc_boarder = 2;
            left_ch_sc_boarder = 0;
            right_ch_sc_boarder = 10;
        }
        Bool_t flag = kFALSE;
        if (
        short_channel_info[channel_number].peak_pos > 11//non-zero signal choosing
        &&  short_channel_info[33].peak_pos > 11
        &&  short_channel_info[32].peak_pos > 11//end

        && abs(short_channel_info[33].peak_pos - short_channel_info[32].peak_pos) < main_sc_boarder// choosing entangled pairs
        && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) > left_ch_sc_boarder // choosing coincidence with NaI counters
        && (short_channel_info[channel_number].peak_pos - short_channel_info[sc_number].peak_pos) < right_ch_sc_boarder) {flag = kTRUE;}
        return flag;
    }

    Bool_t Apply_Amplitude_Saturation_cuts(short_energy_ChannelEntry *short_channel_info, Int_t channel_number, Int_t channel_number_2 = -999)
    {

        Bool_t flag = kFALSE;
        Bool_t flag_ch_2 = kFALSE;
        if (channel_number_2==-999)
            flag_ch_2 = kTRUE;
        else if (short_channel_info[channel_number_2].amp > 200 
                && short_channel_info[channel_number_2].amp < 60000) flag_ch_2 = kTRUE;

        if (short_channel_info[33].amp < 60000
        && short_channel_info[33].amp > 200
        && short_channel_info[32].amp < 60000     
        && short_channel_info[32].amp > 200 
        && short_channel_info[channel_number].amp > 200 
        && short_channel_info[channel_number].amp < 60000
        && short_channel_info[channel_number].integral_in_gate > 1
        && flag_ch_2) flag = kTRUE;
        return flag;
    }


    void double_gauss_fit(TH1F *peak_histo, Float_t &low_cut, Float_t &high_cut, Float_t left_range_width = 1.35, Float_t right_range_width = 1.35, Float_t left_range = 150, Float_t right_range = 350)
    {
            Float_t interval_width = 1.35;

        left_range = peak_histo->GetBinCenter(peak_histo->GetMaximumBin()) - 100;
        right_range = peak_histo->GetBinCenter(peak_histo->GetMaximumBin()) + 100;

        TF1 *gaus_fit = new TF1 ("gaus_func","gaus",left_range,right_range);
        peak_histo->Fit("gaus_func","R");

        TF1 *gaus_2_fit = new TF1 ("gaus_func_2","gaus",
        (gaus_fit->GetParameter(1))-gaus_fit->GetParameter(2),
        (gaus_fit->GetParameter(1))+gaus_fit->GetParameter(2));
        peak_histo->Fit("gaus_func_2","R");

        //mean_int = gaus_2_fit->GetParameter(1);
        //sigma_i = gaus_2_fit->GetParameter(2);
        low_cut = gaus_2_fit->GetParameter(1)-left_range_width*gaus_2_fit->GetParameter(2);
        high_cut = gaus_2_fit->GetParameter(1)+right_range_width*gaus_2_fit->GetParameter(2);
    }
}

#endif

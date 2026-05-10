#include <TH1F.h>
namespace HistPeakToZero
{
    // void GetRanges(TH1F *intitial_hist, Int_t shift_bins, Float_t &xmin, Float_t &xmax, Int_t &nbins)
    // {
    //     nbins = intitial_hist->GetNbinsX();
    //     xmin = (Float_t)(intitial_hist->GetXaxis()->GetXmin() - shift_bins*(intitial_hist->GetXaxis()->GetXmax()-intitial_hist->GetXaxis()->GetXmin())/nbins);
    //     xmax = (Float_t)(intitial_hist->GetXaxis()->GetXmax() - shift_bins*(intitial_hist->GetXaxis()->GetXmax()-intitial_hist->GetXaxis()->GetXmin())/nbins);

    // }
    // void ShiftHist(TH1F *intitial_hist, TH1F *shifted_hist)
    // {
    //     for (Int_t counter = 0; counter < intitial_hist->GetNbinsX()+1; counter++)
    //         shifted_hist->SetBinContent(counter,intitial_hist->GetBinContent(counter));
    // }
    TH1F ShiftedHist(TH1F *intitial_hist, Int_t shift_bins, Int_t sigma = 3, Int_t left_range = 999, Int_t right_range = 999)
    {
        Int_t nbins = intitial_hist->GetNbinsX();

        Float_t xmin = (Float_t)(intitial_hist->GetXaxis()->GetXmin() - shift_bins*(intitial_hist->GetXaxis()->GetXmax()-intitial_hist->GetXaxis()->GetXmin())/nbins);
        Float_t xmax = (Float_t)(intitial_hist->GetXaxis()->GetXmax() - shift_bins*(intitial_hist->GetXaxis()->GetXmax()-intitial_hist->GetXaxis()->GetXmin())/nbins);
        TH1F shifted_hist = TH1F("shifted_hist","",nbins,xmin,xmax);
        for (Int_t counter = 0; counter < intitial_hist->GetNbinsX()+1; counter++)
            shifted_hist.SetBinContent(counter,intitial_hist->GetBinContent(counter));
        if (left_range == 999) left_range = (Int_t)(shifted_hist.GetMean(1)-sigma*shifted_hist.GetStdDev(1));
        if (right_range == 999) right_range = (Int_t)(shifted_hist.GetMean(1)+sigma*shifted_hist.GetStdDev(1));
        shifted_hist.GetXaxis()->SetRangeUser(left_range, right_range);
        cout << shifted_hist.GetMean(1) << " " << shifted_hist.GetStdDev(1) << endl;
        return shifted_hist;
    }
    Int_t find_zero_bin(TH1F *intitial_hist)
    {
        Float_t nbins = intitial_hist->GetNbinsX();
        Float_t xmin = (Float_t)(intitial_hist->GetXaxis()->GetXmin());
        Float_t xmax = (Float_t)(intitial_hist->GetXaxis()->GetXmax());
        Float_t bin_width = (xmax-xmin)/nbins;
        cout << xmax << xmin << bin_width << endl;
        return -((Int_t)(xmin)/bin_width-1);
    }

}
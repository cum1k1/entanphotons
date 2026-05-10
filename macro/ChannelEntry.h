#ifndef CHANNEL_ENTRY_H
#define CHANNEL_ENTRY_H
#include<TTree.h>
#include"constants.h"
//const int MAX_N_SAMPLES = 2048;


struct short_energy_ChannelEntry
{
    Float_t integral_in_gate;
    Float_t peak_pos;
    UShort_t amp; 
    Float_t zl_rms;
	Short_t zero_lvl;


    static TString GetChName(Int_t channel_num)
    {
	    return TString::Format("channel_%i", channel_num);
    }

    TBranch* CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "integral_in_gate/F:peak_pos/F:amp/s:zl_rms/F:zero_lvl/S");
    }

    Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
        return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }

    void Initialize()
    {
        integral_in_gate = 0.;
        peak_pos = 0.;
        amp = 0; 
        zl_rms = 0.;

    }

    
};

struct diff_short_energy_ChannelEntry
{
    Short_t min_diff;
    Short_t min_diff_time;
    Short_t max_diff;
    Short_t max_diff_time;

    static TString GetChName(Int_t channel_num)
    {
	    return TString::Format("diff_channel_%i", channel_num);
    }

    TBranch* CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "min_diff/S:min_diff_time/S:max_diff/S:max_diff_time/S");
    }

    Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
        return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }

    void Initialize()
    {
        min_diff = 0;
        min_diff_time = 0;
        max_diff = 0;
        max_diff_time = 0;
    }    
    
};

struct ChannelEntry {

    Short_t wf_size;
    Short_t wf[MAX_N_SAMPLES];
    

    static TString GetChName(Int_t channel_num)
    {
	return TString::Format("channel_%i", channel_num);
    }

    Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
	return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }
    
    void Initialize()
    {
        for (int i = 0; i < sizeof(wf)/sizeof(wf[0]); i++) wf[i] = 0;
        wf_size = 0;
    }

    void SplineWf()
    {
        Float_t wf1[MAX_N_SAMPLES] = {0};
        const Int_t SplineWidth = 2;
        for (Int_t i = 0; i < wf_size; i++)
        {
            Int_t il=i-SplineWidth; Int_t ir=i+SplineWidth;
            if (il<0) il=0;
            if (ir>wf_size-1) ir=wf_size-1;
            Float_t counter = 0;
            for (Int_t in = il; in <=ir; in++) {wf1[i]+=wf[in];counter++;}
            wf1[i]/=counter;
        }
        for (Int_t i = 0; i < wf_size; i++) wf[i] = wf1[i];
    }

    void DiffWf()
    {
        const Float_t Diff_window = 4;
        Short_t wf1[MAX_N_SAMPLES] = {0};
        for (Int_t i = 0; i < wf_size; i++)
        {
            Int_t il=i-Diff_window; Int_t ir=i+Diff_window;
            if (il<0) il=0;
            if (ir>wf_size-1) ir=wf_size-1;
            wf1[i]=(Short_t)((wf[ir]-wf[il])/(ir-il));
        }
        for (Int_t i = 0; i < wf_size; i++) wf[i] = wf1[i];
    }    

    void FindDiffWfPars(Short_t &min_diff, Short_t &min_time, Short_t &max_diff, Short_t &max_time, Int_t GATE_BEG = 55, Int_t GATE_END = 76)
    {
        //Int_t vmin = numeric_limits<Int_t>::max();
        //Int_t vmax = numeric_limits<Int_t>::min();

        for (Short_t s=GATE_BEG; s < GATE_END; ++s) {
            Short_t v = wf[s];
            if (v < min_diff) {
                min_diff = v;
                min_time = 16*s;
            }
            if (v > max_diff) {
                max_diff = v;
                max_time = 16*s;
            }      
        }
    }

    Int_t Get_Zero_Level(Int_t GATE_BEG)
    {
        const Int_t interv_num = 1;
        int zero_lvl = 0;
        int best_spread = -1;
        for (int i=0; i < interv_num; ++i) {
            int vmin = numeric_limits<int>::max();
            int vmax = numeric_limits<int>::min();
            int sum = 0;
            for (int s=GATE_BEG/interv_num * i; s < GATE_BEG/interv_num * (i+1); ++s) {
                int v = wf[s]; 
                sum += v;
                if (v < vmin) vmin = v;
                if (v > vmax) vmax = v;
            }
            int spread = vmax - vmin;
            if (best_spread < 0) best_spread = spread;
            if (spread <= best_spread) {
                best_spread = spread;
                zero_lvl = sum / (GATE_BEG/interv_num);
            }
        }
        return zero_lvl;
    }
    Float_t Get_Zero_Level_RMS(Int_t GATE_BEG)
    {
        const Int_t interv_num = 1;
        int zero_lvl = 0;
        Float_t best_spread = -1;
        Float_t rms_zl = -1;
        for (Int_t i=0; i < interv_num; ++i) {
            Int_t vmin = numeric_limits<int>::max();
            Int_t vmax = numeric_limits<int>::min();
            Float_t sum = 0; Float_t sumsquare = 0; Float_t sum_counter = 0;
            for (Int_t s=GATE_BEG/interv_num * i; s < GATE_BEG/interv_num * (i+1); ++s) {
                Int_t v = wf[s]; 
                sum += (Float_t)v;
                //sumsquare += v*v;
                sum_counter++;
            }
            sum /=sum_counter;
            sumsquare = 0.;
            
            for (Int_t s=GATE_BEG/interv_num * i; s < GATE_BEG/interv_num * (i+1); ++s) {
                sumsquare += (Float_t)(wf[s] - sum)*(wf[s] - sum)/sum_counter;
            }            
            rms_zl = sqrt(sumsquare);
            //rms_zl = (Float_t)(sumsquare/sum_counter-pow((sum/sum_counter),2));

            if (best_spread < 0) best_spread = rms_zl;
            //cout << best_spread<<endl<<"rf3rv"<<endl;
            // if (rms_zl <= best_spread) {
            //     best_spread = rms_zl;
            //     zero_lvl = sum / sum_counter;
            // }
        }
        return best_spread;
    }

    Float_t Get_Integral_in_Gate(Int_t zero_level, Int_t GATE_BEG, Int_t GATE_END, Int_t pp =200)
    {
        Float_t gateInteg = 0;
	//Int_t GATE_BEG1 = pp -40;
	//Int_t GATE_END1 = pp+60;
	//if (GATE_BEG1 < 0) GATE_BEG1=0;
	//if (GATE_END1 > GATE_END) GATE_END1 = GATE_END;
        for (int s=GATE_BEG; s < GATE_END+1; ++s) {
            	if ((float)zero_level > (float)wf[s]) 
			gateInteg +=  (float)zero_level - (float)wf[s] ;
		//if ((float)zero_level < (float)wf[s] && s > pp) break;


        }
        return gateInteg;
    }
    
    Short_t Get_peak_position(Int_t zero_lvl)
    {
        Int_t amp = 0;
        Short_t peakPos = 0;
        for (int s=0; s < wf_size; ++s) {
            int v = wf[s] - zero_lvl;
            if (v < amp) {
                amp = v;
                peakPos = s;
            }
        }
        return peakPos;
    }
    Float_t Get_peak_position_gauss(Int_t zero_level, Int_t inv_amp, Int_t CH_GATE_END = -1000)
    {
        if ( CH_GATE_END ==-1000) CH_GATE_END = wf_size;
        if (wf_size < CH_GATE_END) CH_GATE_END = wf_size;
        Float_t peak_search = 0.;
        Float_t ampl_sum = 0;
        for (Int_t s= 0; s < CH_GATE_END; ++s) {
            Int_t v = zero_level - wf[s];
            if (v > inv_amp*0.1)
            {
                ampl_sum += (Float_t) v;
                peak_search+= (Float_t) v*s;
            }
        }
        peak_search /= ampl_sum;
    return 16.0*peak_search;
    }
////////////
    UShort_t Get_Amplitude(Int_t zero_level,Int_t GATE_BEG = 0, Int_t CH_GATE_END = -1000)
    {
        if ( CH_GATE_END ==-1000) CH_GATE_END = wf_size;
        if (wf_size == 0) {GATE_BEG = 0; CH_GATE_END = 0;}

        UShort_t amp = numeric_limits<UShort_t>::min();
        for (int s=GATE_BEG; s < CH_GATE_END; ++s) {
            Int_t v =  (Int_t)zero_level - (Int_t)wf[s];
            if (v > amp) {amp = v;}
        }
        return amp;
    }

    // void Calculate_waveform(short_energy_ChannelEntry &result_event, Int_t ch_num, Int_t gate_beg, Int_t gate_end, Int_t gate_maximum_beg, Int_t gate_maximum_end,  Bool_t IsFIT = true, Double_t first_fit_harmonic = 0., Double_t second_fit_harmonic = 0., TString FIT_QA_mode = "", TObjArray *check_fit_arr = NULL, Float_t *fitQA_arg_arr = NULL);
    
    // void Calculate_fit_harmonics(fit_harmonics_struct &result_fit_event, Int_t &event_counter, Int_t event_num, Int_t ch_num, Int_t gate_beg, Int_t gate_end, Int_t gate_maximum_beg, Int_t gate_maximum_end,TString source_path, TString run_name, TObjArray *check_fit_arr = NULL);

};

//##################
//#################


// struct short_ChannelEntry
// {
//     Int_t integral_in_gate;
//     Short_t peak_pos;
//     UShort_t amp; 

//     static TString GetChName(Int_t channel_num)
//     {
// 	return TString::Format("channel_%i", channel_num);
//     }

//     TBranch* CreateBranch(TTree *tree, Int_t channel_num)
//     {
// 	    return tree->Branch(GetChName(channel_num).Data(), this, "integral_in_gate/I:peak_pos/S:amp/s");
//     }
//         Int_t SetBranch(TTree *tree, Int_t channel_num)
//     {
// 	return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
//     }
// };


struct mini_tree_nrg
{
    Float_t EdepIntermediate;
    Float_t EdepScat0;
    Float_t EdepScat1;
    Float_t EdepDet0;
    Float_t EdepDet1;  
    Short_t DetNum0;
    Short_t DetNum1; 


    static TString GetBrName()
    {
	    return TString::Format("MiniTree");
    }

    TBranch* CreateBranches(TTree *tree)
    {
	    return tree->Branch(GetBrName().Data(), this, "EdepIntermediate/F:EdepScat0/F:EdepScat1/F:EdepDet0/F:EdepDet1/F:DetNum0/S:DetNum1/S");
    }

    Int_t SetBranch(TTree *tree)
    {
        return tree->SetBranchAddress(GetBrName().Data(), this);
    }

    Int_t Initialize()
    {
        EdepIntermediate = 0;
        EdepScat0 = 0;
        EdepScat1 = 0;
        EdepDet0 = 0;
        EdepDet1 = 0;  
        DetNum0 = 0;
        DetNum1 = 0; 
        return 1;
    }
};

struct mini_tree_time
{
    Float_t TimeIntermediate;
    Float_t TimeScat0;
    Float_t TimeScat1;
    Float_t TimeDet0;
    Float_t TimeDet1;

    static TString GetBrName()
    {
	    return TString::Format("TimeTree");
    }

    TBranch* CreateBranches(TTree *tree)
    {
	    return tree->Branch(GetBrName().Data(), this, "TimeIntermediate/F:TimeScat0/F:TimeScat1/F:TimeDet0/F:TimeDet1/F");
    }

    Int_t SetBranch(TTree *tree)
    {
        return tree->SetBranchAddress(GetBrName().Data(), this);
    }

    Int_t Initialize()
    {
        TimeScat0 = 0;
        TimeScat1 = 0;
        TimeDet0 = 0;
        TimeDet1 = 0;
        TimeIntermediate = 0;
        return 1;
    }
};

struct mini_tree
{

    Float_t EdepScat0;
    Float_t EdepScat1;
    Float_t EdepDet0;
    Float_t EdepDet1;  
    Float_t EdepWeak;
    Short_t DetNum0;
    Short_t DetNum1; 
    Float_t TimeScat0;
    Float_t TimeScat1;
    Float_t TimeDet0;
    Float_t TimeDet1;
    Float_t TimeWeak;
    std::map<TString, Float_t*> branchfloat = 
    {
        {"EdepScat0", &EdepScat0},
        {"EdepScat1", &EdepScat1},
        {"EdepDet0", &EdepDet0},
        {"EdepDet1", &EdepDet1},     
        {"EdepWeak", &EdepWeak}       
    }; 
    std::map<TString, Short_t*> branchshort = 
    {
        {"DetNum0", &DetNum0},
        {"DetNum1", &DetNum1}
    }; 
    std::map<TString, Float_t*> branchfloat2 = 
    {
        {"TimeScat0", &TimeScat0}
        ,{"TimeScat1", &TimeScat1}
        ,{"TimeDet0", &TimeDet0}
        ,{"TimeDet1", &TimeDet1}
        ,{"TimeWeak", &TimeWeak}
    };

    Int_t CreateBranches(TTree* mini_tree)
    {
        for (auto& val:branchfloat)
            mini_tree->Branch(val.first,val.second,val.first+"/F");
        for (auto& val:branchshort)
            mini_tree->Branch(val.first,val.second,val.first+"/S");        
        for (auto& val:branchfloat2)
            mini_tree->Branch(val.first,val.second,val.first+"/F");
        return 1;
    }

    void SetBranch(TTree* mini_tree)
    {

        for (auto& val:branchfloat)
            mini_tree->SetBranchAddress(val.first,val.second);
        for (auto& val:branchshort)
            mini_tree->SetBranchAddress(val.first,val.second);        
        for (auto& val:branchfloat2)
            mini_tree->SetBranchAddress(val.first,val.second);

    }

    Int_t Initialize()
    {
        for (auto& val:branchfloat)        
        {
            val.second = 0;
        }
        for (auto& val:branchshort)
        {
            val.second = 0;
        }
        return 1;
    }

};
#endif CHANNEL_ENTRY_H

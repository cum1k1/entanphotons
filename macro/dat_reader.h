
#include<TTree.h>
#include<TString.h>
#include<TTree.h>
const int MAX_N_SAMPLES = 2048;

struct ecal_data
{
    Int_t wf[MAX_N_SAMPLES];
    Short_t wf_size;

    void Initialize()
    {
        for (Int_t i = 0; i < MAX_N_SAMPLES; i++) wf[i] = 0;
        wf_size = 60;
    }

    void ReadWf(string data_line)
    {
        Int_t counter = 0;
        stringstream linestream(data_line);
        Int_t value;
        while (linestream >> value)
        {
            wf[counter] = value;
            //cout << value << " ";
            counter++;
        }

    }
    void OnlyWf(Int_t length)
    {
        for (Int_t i = 0; i < length-6; i++) wf[i] = wf[i+6];
    }
    void ShiftWf(Int_t length,Int_t zl)
    {
        for (Int_t i = 0; i < length-6; i++) wf[i] = wf[i]-zl;
    }

    
    static TString GetChName(Int_t channel_num)
    {
	return TString::Format("channel_%i", channel_num);
    }

    Int_t SetBranch(TTree *tree, Int_t channel_num)
    {
	return tree->SetBranchAddress(GetChName(channel_num).Data(), this);
    }

    void InvertSignal()
    {
        for (Int_t i = 0; i < wf_size; i++)
        {
            wf[i] = -wf[i];
        }
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

    void FindDiffWfPars(Short_t &min_diff, Short_t &min_time, Short_t &max_diff, Short_t &max_time, Int_t GATE_BEG = 15, Int_t GATE_END = 36)
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
    Int_t Get_Zero_Level(Int_t GATE_BEG=20)
    {
        const Int_t interv_num = 1;
        int zero_lvl = 0;
        int best_spread = -1;
        float sum = 0;
        for (int i=0; i < GATE_BEG; ++i) {
	sum += wf[i];

	}
      	zero_lvl = (int)(sum / (GATE_BEG));            
        return zero_lvl;
    }
    Float_t Get_Zero_Level_RMS(Int_t GATE_BEG=20)
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
                sum_counter++;
            }
            sum /=sum_counter;
            sumsquare = 0.;
            
            for (Int_t s=GATE_BEG/interv_num * i; s < GATE_BEG/interv_num * (i+1); ++s) {
                sumsquare += (Float_t)(wf[s] - sum)*(wf[s] - sum)/sum_counter;
            }            
            rms_zl = sqrt(sumsquare);
            if (best_spread < 0) best_spread = rms_zl;

        }
        return best_spread;
    }

    Float_t Get_Integral_in_Gate(Int_t zero_level, Int_t GATE_BEG =20, Int_t GATE_END=36)
    {
        Float_t gateInteg = 0;
        for (int s=GATE_BEG; s < GATE_END+1; ++s) {
            if (wf[s] < zero_level) gateInteg +=  (float)(zero_level - wf[s]) ;
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
    return peak_search;
    }
////////////
    Int_t Get_Amplitude(Int_t zero_level,Int_t GATE_BEG = 0, Int_t CH_GATE_END = -1000)
    {
        if ( CH_GATE_END ==-1000) CH_GATE_END = wf_size;
        if (wf_size == 0) {GATE_BEG = 0; CH_GATE_END = 0;}

        Int_t amp = numeric_limits<UShort_t>::min();
        for (int s=GATE_BEG; s < CH_GATE_END; ++s) {
            Short_t v =  zero_level - wf[s];
            if (v > amp) {amp = v;}
        }
        return amp;
    }
    
};


struct short_energy_ChannelEntry
{
    Float_t integral_in_gate;
    Float_t peak_pos;
    UShort_t amp; 
    Float_t zl_rms;

    static TString GetChName(Int_t channel_num)
    {
	    return TString::Format("channel_%i", channel_num);
    }
    TBranch* CreateBranch(TTree *tree, Int_t channel_num)
    {
	    return tree->Branch(GetChName(channel_num).Data(), this, "integral_in_gate/F:peak_pos/F:amp/s:zl_rms");
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

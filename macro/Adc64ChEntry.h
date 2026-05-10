#ifndef ADC64CHENTRY_H
#define ADC64CHENTRY_H

#include <TRefArray.h>

const int MAX_N_SAMPLES = 2048;

class Adc64ChEntry : public TObject {
public:
    Int_t integral_in_gate;
    Short_t zero_level;
    Short_t peak_pos;
    Short_t amp;
    // 16-bit gap
    Int_t wf_size;
    Short_t wf_plug[MAX_N_SAMPLES]; //!
    Short_t *wf;                    //[wf_size]

    Adc64ChEntry() : wf(0) {}
    Adc64ChEntry(void *dummy) : wf(wf_plug) {} 

ClassDef(Adc64ChEntry, 1)
};


const char *EAST_DETECTORS[] = {
/* 1:   */    "channel_0",
/* 2:   */    "channel_1",
/* 3:   */    "channel_2",
/* 4:   */    "channel_3",
/* 5:   */    "channel_4",
/* 6:   */    "channel_5",
/* 7:   */    "channel_6",
/* 8:   */    "channel_7",
/* 9:   */    "channel_8",
/* 10:  */    "channel_9",
/* 11:  */    "channel_10",
/* 12:  */    "channel_11",
/* 13:  */    "channel_12",
/* 14:  */    "channel_13",
/* 15:  */    "channel_14",
/* 16:  */    "channel_15",
}; 

class TTree;

class ArmEntry : public TObject {
public:
    TRefArray detectors;

    ArmEntry(TTree *adc64);

ClassDef(ArmEntry, 1)
};

#endif
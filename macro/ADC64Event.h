#ifndef ADC64_EVENT_H
#define ADC64_EVENT_H

#include <Rtypes.h>
#include <vector>

struct ADC64Frame 
{
    Short_t     wf_size;
    Short_t*    wf;         //[wf_size]
    bool        isOwner;    //!

    ADC64Frame(bool=false);
    ~ADC64Frame();

ClassDef(ADC64Frame, 1)
};

struct ADC64Event
{
    std::vector<ADC64Frame*> frames;
    Long64_t    ts;

    ADC64Event(int=0);
    ~ADC64Event();

    void reset();
    ADC64Frame* operator[](int ch) { return frames[ch]; }

ClassDef(ADC64Event, 1)
};

#endif
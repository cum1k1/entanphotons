#include "ADC64Event.h"
#include <chrono>

using namespace std::chrono;

const int MAX_N_SAMPLES = 2048;

ADC64Frame::ADC64Frame(bool doInit)
    : wf_size(0),
    wf(doInit ? new Short_t[MAX_N_SAMPLES] : nullptr),
    isOwner(doInit)
{}

ADC64Frame::~ADC64Frame()
{
    if (isOwner) delete [] wf;
}


ADC64Event::ADC64Event(int n_channels)
{
    for (int ch=0; ch < n_channels; ++ch) frames.push_back( new ADC64Frame(true) );
}

ADC64Event::~ADC64Event()
{
    for (auto fr : frames) delete fr;
}

void ADC64Event::reset()
{
    for (auto fr : frames) fr->wf_size = 0;
    ts = duration_cast<seconds>(system_clock::now().time_since_epoch()).count();
}
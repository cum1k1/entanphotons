#ifndef CHSH_CALCULATOR_H
#define CHSH_CALCULATOR_H

namespace CHSH
{



    void count_coincidences(Int_t n, Int_t NumEvents[16][16], Int_t &N_parallel, Int_t &N_perpendicular, Int_t &N_mixed_1, Int_t &N_mixed_2, TString z = "clockwise")
    {
            N_parallel=0;
            N_perpendicular = 0;
            N_mixed_1 = 0;
            N_mixed_2 = 0;  

        for (Int_t a = 0; a < 16; a++)
        {
            Int_t b = a+n+48;
            Int_t a_1 = a+4;
            Int_t b_1 = b+4;
            while (a > 15) a=a-16;
            while (b > 15) b=b-16;
            while (a_1> 15) a_1=a_1-16;
            while (b_1 > 15) b_1=b_1-16;    
            N_parallel += NumEvents[a][b];
            N_perpendicular += NumEvents[a_1][b_1];
            N_mixed_1 += NumEvents[a][b_1];
            N_mixed_2 += NumEvents[a_1][b];        
        }
        if (z=="both")
            for (Int_t a = 0; a < 16; a++)
            {
                Int_t b = a-n+48;
                Int_t a_1 = a+4;
                Int_t b_1 = b+4;
                while (a > 15) a=a-16;
                while (b > 15) b=b-16;
                while (a_1> 15) a_1=a_1-16;
                while (b_1 > 15) b_1=b_1-16;    
                N_parallel += NumEvents[a][b];
                N_perpendicular += NumEvents[a_1][b_1];
                N_mixed_1 += NumEvents[a][b_1];
                N_mixed_2 += NumEvents[a_1][b];        
            }
    }

    Float_t global_E_coeff(Int_t n, Int_t NumEvents[16][16], TString z = "clockwise")
    {
        Int_t N_parallel = 0; Int_t N_perpendicular = 0;Int_t N_mixed_1 = 0;Int_t N_mixed_2 = 0;
        count_coincidences(n,NumEvents, N_parallel, N_perpendicular, N_mixed_1, N_mixed_2,z);
        return (float)(N_parallel+N_perpendicular-N_mixed_1-N_mixed_2)/(N_parallel+N_perpendicular+N_mixed_1+N_mixed_2)
        ;
    }
    Float_t global_sqr_E_error(Int_t n,Int_t NumEvents[16][16], TString z = "clockwise") //calculating square error of correlation coefficient
    {
        Int_t N_parallel = 0; Int_t N_perpendicular = 0;Int_t N_mixed_1 = 0;Int_t N_mixed_2 = 0;

        count_coincidences(n,NumEvents, N_parallel, N_perpendicular, N_mixed_1, N_mixed_2,z);
        return (float)4.0/pow((N_parallel+N_perpendicular+N_mixed_1+N_mixed_2),4)*(pow(N_mixed_1+N_mixed_2,2)*(N_parallel+N_perpendicular)+pow(N_parallel+N_perpendicular,2)*(N_mixed_1+N_mixed_2));
    }
    Float_t global_calculate_CHSH(Int_t n,Int_t NumEvents[16][16], TString z = "clockwise")
    {
        if (z == "counterclockwise") n = -n;
        return
        3*global_E_coeff(n, NumEvents,z) - global_E_coeff(3*n,NumEvents,z)
        ;
    }

    Float_t global_calculate_CHSH_error(Int_t n,Int_t NumEvents[16][16], TString z = "clockwise")
    {
        if (z == "counterclockwise") n = -n;

        return
        global_sqr_E_error(n, NumEvents,z)*3 + global_sqr_E_error(3*n, NumEvents,z)
        ;
    }

/////////////////////////////////////////
/////////////////////////////////////////local_CHSH_for_average
    void local_count_coincidences(Int_t a, Int_t b, Int_t NumEvents[16][16], Int_t &N_parallel, Int_t &N_perpendicular, Int_t &N_mixed_1, Int_t &N_mixed_2, TString z = "clockwise")
    {
        a += 48;
        b += 48;
        Int_t a_1 = a+4;
        Int_t b_1 = b+4;
        while (a > 15) a=a-16;
        while (b > 15) b=b-16;
        while (a_1> 15) a_1=a_1-16;
        while (b_1 > 15) b_1=b_1-16;    
        N_parallel=NumEvents[a][b];
        N_perpendicular = NumEvents[a_1][b_1];
        N_mixed_1 = NumEvents[a][b_1];
        N_mixed_2 = NumEvents[a_1][b];
    }

    Float_t E_coeff(Int_t NumEvents[16][16], Int_t a, Int_t b, TString z = "clockwise")
    {
        Int_t N_parallel = 0;
        Int_t N_perpendicular = 0;
        Int_t N_mixed_1 = 0;
        Int_t N_mixed_2 = 0;
        local_count_coincidences(a,b,NumEvents, N_parallel, N_perpendicular, N_mixed_1, N_mixed_2);
        return (float)(N_parallel+N_perpendicular-N_mixed_1-N_mixed_2)/(N_parallel+N_perpendicular+N_mixed_1+N_mixed_2)
        ;
    }

    Float_t sqr_E_error(Int_t NumEvents[16][16], Int_t a, Int_t b, TString z = "clockwise") //calculating square error of correlation coefficient
    {
        Int_t N_parallel = 0;
        Int_t N_perpendicular = 0;
        Int_t N_mixed_1 = 0;
        Int_t N_mixed_2 = 0;
        local_count_coincidences(a,b,NumEvents, N_parallel, N_perpendicular, N_mixed_1, N_mixed_2);
        return 
        (float)4.0/pow((N_parallel+N_perpendicular+N_mixed_1+N_mixed_2),4)*(pow(N_mixed_1+N_mixed_2,2)*(N_parallel+N_perpendicular)+pow(N_parallel+N_perpendicular,2)*(N_mixed_1+N_mixed_2))
    ;
    }

    Float_t calculate_local_CHSH(Int_t NumEvents[16][16], Int_t a_angle, Int_t angle_between_counters, TString z = "clockwise")
    {
        if (z == "counterclockwise") angle_between_counters = -angle_between_counters;
        return
        E_coeff(NumEvents,a_angle,a_angle+angle_between_counters,z) -
        E_coeff(NumEvents,a_angle,a_angle+3*angle_between_counters,z) +
        E_coeff(NumEvents,a_angle+2*angle_between_counters,a_angle+angle_between_counters,z) +
        E_coeff(NumEvents,a_angle+2*angle_between_counters,a_angle+3*angle_between_counters,z);
    }

    Float_t calculate_local_CHSH_error(Int_t NumEvents[16][16], Int_t a_angle, Int_t angle_between_counters, TString z = "clockwise")
    {
        if (z == "counterclockwise") angle_between_counters = -angle_between_counters;
        return
        sqr_E_error(NumEvents,a_angle,a_angle+angle_between_counters,z) +
        sqr_E_error(NumEvents,a_angle,a_angle+3*angle_between_counters,z) +
        sqr_E_error(NumEvents,a_angle+2*angle_between_counters,a_angle+angle_between_counters,z) +
        sqr_E_error(NumEvents,a_angle+2*angle_between_counters,a_angle+3*angle_between_counters,z);
    }

    Int_t true_number(Int_t a)
    {
        a+=48;
        while(a > 15) a = a-16;
        return a;
    }
}
#endif
    
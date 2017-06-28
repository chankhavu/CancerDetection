#ifndef OTSU_THRESHOLD_HEADER
#define OTSU_THRESHOLD_HEADER


void CreateHistogram(int** GrayImg, int* Hist, int Width, int Height, int* MinValue, int* MaxValue){
    *MinValue = GrayImg[0][0];
    *MaxValue = GrayImg[0][0];

    for (int i=0; i<256; i++)
        Hist[i] = 0;

    for (int i=0; i<Width; i++)
        for (int j=0; j<Height; j++){
            Hist[ GrayImg[i][j] ]++;
            if (GrayImg[i][j] > *MaxValue)
                *MaxValue = GrayImg[i][j];
            if (GrayImg[i][j] < *MinValue)
                *MinValue = GrayImg[i][j];
        }
}

int OtsuThreshold(int** GrayImg, int Width, int Height){
    // First, create histogram
    int Histogram[256], MinValue, MaxValue;
    CreateHistogram(GrayImg, Histogram, Width, Height, &MinValue, &MaxValue);

    int m = 0, n = 0;
    for (int t = MinValue; t <= MaxValue; t++)
        m += (t - MinValue) * Histogram[t], n += Histogram[t];

    float MaxSigma = -1.0;
    int Threshold = 0;
    int Alpha1 = 0, Beta1 = 0;

    for (int t = MinValue; t <= MaxValue; t++){
        Alpha1 += (t - MinValue) * Histogram[t];
        Beta1 += Histogram[t];

        // possibility of 1st class, w2 = 1-w1
        float w1 = (float)Beta1 / (1. * n);

        // a = a1 - a2, where a1 and a2 are arithmetic mean for classes 1 and 2
        float a = (float)Alpha1 / (1. * Beta1) - (float)(m - Alpha1) / (1. * (n - Beta1));

        // now calculate sigma
        float sigma = w1 * (1.0 - w1) * a * a;

        if (sigma > MaxSigma){
            MaxSigma = sigma;
            Threshold = t;
        }
    }

    return Threshold;
}


#endif

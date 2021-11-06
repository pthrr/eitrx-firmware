#include <stdint.h>
#include "fft.h"

#define SWAP(a,b) tempr=(a); (a)=(b); (b)=tempr

const double kWTEMP_TABLE[9] = {
  #include "wtemp_512.inc"
};

const double kWPR_TABLE[9] = {
  #include "wpr_512.inc"
};

const double kWPI_TABLE[9] = {
  #include "wpi_512.inc"
};

void diffft2(float* data,
             const uint16_t kN,
             const int8_t kISIGN)
{
  float wtemp;
  float wr;
  float wi;
  float wpr;
  float wpi;
  float tempr;
  float tempi;
  uint8_t k = 0;
  uint16_t i = 0;
  uint16_t j = 0;
  //uint16_t n = 0;
  uint16_t m = 0;
  uint16_t istep;
  uint16_t mmax;

  data = data - 1; // prepare for iterations
  //n = kN;
  mmax = kN/2;

  // calculate the FFT
  while (mmax >= 2)
  {
      istep = mmax << 1;
      // TODO switch tables for kISIGN
      wtemp = kWTEMP_TABLE[k];
      wpr = kWPR_TABLE[k];
      wpi = kWPI_TABLE[k];
      k += 1;
      wr = 1.0;
      wi = 0.0;

      for (m = 1; m < mmax; m += 2)
      {
          for (i = m; i <= kN; i += istep)
          {
              j = i + mmax;
              tempr = data[i];
              tempi = data[i+1];
              data[i] = data[i] + data[j];
              data[i+1] = data[i+1] + data[j+1];
              data[j] = tempr - data[j];
              data[j+1] = tempi - data[j+1];
              tempr = wr*data[j] - wi*data[j+1];
              tempi = wr*data[j+1] + wi*data[j];
              data[j] = tempr;
              data[j+1] = tempi;
          }

          wtemp = wr;
          wr += wtemp*wpr - wi*wpi;
          wi += wtemp*wpi + wi*wpr;
      }

      mmax = mmax/2;
  }

  // do the bit-reversal
  j = 1;

  for (i = 1; i < kN; i += 2)
  {
      if (j > i)
      {
          SWAP(data[j], data[i]);
          SWAP(data[j+1], data[i+1]);
      }

      m = kN >> 1;

      while (m >= 2 && j > m)
      {
          j -= m;
          m >>= 1;
      }

      j += m;
  }
} // end of diftt2()

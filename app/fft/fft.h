#ifndef _FFT_H
#define _FFT_H

#include <stdint.h>

void diffft2(float* data,
             const uint16_t kN,
             const int8_t kISIGN);

#endif // _FFT_H

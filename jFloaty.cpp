// jFloaty.cpp : Decode a JPEG direct to RGBFFF(float32)
//
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#define STBI_NO_SIMD
#define STB_IMAGE_IMPLEMENTATION
#define STBI_ONLY_JPEG
#include "stb_image.h"

static int desired_channels = 1;

static int fltVecCmp(const void* l_, const void* r_) {
  const float *l = (const float*)l_;
  const float *r = (const float*)r_;
  for (unsigned i = 0; i < desired_channels; ++i) {
    if (r[i] > l[i]) return 1;
    if (l[i] > r[i]) return-1;
  }
  return 0;
}


int main()
{
    std::cout << "Hello World!\n";

#if 0
    for (unsigned i = 0; i < 10; ++i) {
        short idct_input[64];
        for (unsigned j = 0; j < 64; ++j) {
            idct_input[j] = rand() % 64 - 32;
        }

        stbi_uc out_fix[64];
        stbi__idct_block(out_fix, 8, idct_input);

        float out_f32[64];
        stbi__idct_block_f32(out_f32, 8, idct_input);

        float maxDiff = 0;
        for (unsigned j = 0; j < 64; ++j) {
            float diff = fabsf(out_fix[j] / 255.0f - out_f32[j]);
            if (diff > maxDiff) {
                maxDiff = diff;
            }
        }

        std::cout << "Run:" << i << " maxDiff:" << maxDiff << std::endl;
    }
#endif

    // --- Reading from a file (mode: "r" for read) ---
    FILE* fp = fopen("input.jpg", "rb");
    if (fp == NULL) {
        printf("Error opening file for reading!\n");
        exit(__LINE__);
    }

    int w, h;
    int channels_in_file;

    std::cout << "Checking channel count:";
    if (!stbi_info_from_file(fp, &w, &h, &channels_in_file)) {
      std::cerr << "stbi_info_from_file() failure" << std::endl;
      exit(__LINE__);
    }
    std::cout << channels_in_file << std::endl;

    std::cout << "Reading data from file...";
    desired_channels = channels_in_file;
    float* pixels = stbi_loadf_from_file(fp, &w, &h, &channels_in_file, desired_channels);
    if (!pixels) {
      std::cerr << "stbi_loadf_from_file() failure" << std::endl;
      exit(__LINE__);
    }
    fclose(fp); // Close the file
    std::cout << "done" << std::endl;

    if (pixels) {
      std::cout << "width:" << w << std::endl;
      std::cout << "height:" << h << std::endl;

      fp = fopen("output.fff.data", "wb");
      assert(fp);
      fwrite(pixels, sizeof(float), w * h * desired_channels, fp);
      fclose(fp);
    }
    else {
        std::cout << stbi_failure_reason() << std::endl;
        return -1;
    }

    std::cout << "Unique pixels:";
    qsort(pixels, w * h, sizeof(float) * desired_channels, fltVecCmp);
    unsigned pcount = 1;
    for (unsigned i = 1; i < w * h; ++i) {
      const float* pixel = pixels + i * desired_channels;
      pcount += !memcmp(pixel, pixel - desired_channels, sizeof(float) * desired_channels);
    }
    std::cout << pcount << std::endl;

    return 0;
}
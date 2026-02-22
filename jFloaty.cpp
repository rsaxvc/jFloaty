// jFloaty.cpp : Decode a JPEG direct to RGBFFF(float32)
//
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#define STBI_NO_SIMD
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

static int desired_channels = 1;

static int fltVecCmp(const void* l_, const void* r_) {
  const float* l = (const float*)l_;
  const float* r = (const float*)r_;
  for (int i = 0; i < desired_channels; ++i) {
    if (r[i] > l[i]) return 1;
    if (l[i] > r[i]) return-1;
  }
  return 0;
}

static void die(int line) {
  std::cerr << "die() called at line #" << line << std::endl;
  exit(line);
}

static void idct_test() {
  std::cout << "running inverse DCT comparison" << std::endl;
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

    std::cout << "\tRun:" << i << " maxDiff:" << maxDiff << std::endl;
  }
}

static void usage(const char* progName) {
  std::cerr << progName << ": IEEE 32-bit Floating-Point JPEG codec" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Mental Model:" << std::endl;
  std::cerr << "\tArguments are processed in left to right order." << std::endl;
  std::cerr << "\tThe program works with one image-buffer at a time." << std::endl;
  std::cerr << "\tLoading populates the image-buffer." << std::endl;
  std::cerr << "\tStoring writes the image-buffer to a file." << std::endl;
  std::cerr << std::endl;
  std::cerr << "Error handling: immediate, and fatal" << std::endl;
  std::cerr << std::endl;
  std::cerr << "Arguments:" << std::endl;
  std::cerr << "\t-i <filename> #load file into buffer" << std::endl;
  std::cerr << "\t-s #generate some statistics" << std::endl;
  std::cerr << "\t-o <filename> #safe buffer to file" << std::endl;
  std::cerr << "\t-o16 <filename> #safe buffer to 16-bit file" << std::endl;
  std::cerr << std::endl;
  std::cerr << "\t-w <width> #Set image width for raw-float loading" << std::endl;
  std::cerr << std::endl;
  std::cerr << "\t-idct #run idct self test" << std::endl;
  std::cerr << std::endl;
  std::cerr << "\t-q 0-99 #Set JPEG encoder quality factor" << std::endl;
  std::cerr << std::endl;
  std::cerr << std::endl;
  std::cerr << "Example: RGBFFF -> JPG" << std::endl;
  std::cerr << "\t" << progName << " -w 1280 -i input.fff.data -q 95 -o output.jpg" << std::endl;
  std::cerr << "Example: JPG -> RGBFFF" << std::endl;
  std::cerr << "\t" << progName << " -i input.jpg -o output.fff.data" << std::endl;
  std::cerr << "Example: JPG -> PNG" << std::endl;
  std::cerr << "\t" << progName << " -i input.jpg -o output.png" << std::endl;
  std::cerr << std::endl;
  std::cerr << std::endl;
}

static bool strEndsWith(const char* haystack, const char* needle) {
  const auto haystackLen = strlen(haystack);
  const auto needleLen = strlen(needle);
  if (haystackLen < needleLen) {
    return false;
  }
  return 0 == strcmp(haystack + haystackLen - needleLen, needle);
}

static float * rot90(float * in, int & w, int & h, int c) {
  float* out = (float*)malloc(sizeof(float) * w * h * c);
  if (!out) {
    std::cerr << "out of memory" << std::endl;
    die(__LINE__);
  }

  for (int x = 0; x < w; ++x) {
    for (int y = 0; y < h; ++y) {
      memcpy(out + c * ((w-x-1) * h + y), in + c * (y * w + x), sizeof(*in) * c);
    }
  }
  std::swap(w, h);
  free(in);
  return out;
}

static float* rot180(float* in, int & w, int & h, int c) {
  return rot90(rot90(in, w, h, c), w, h, c);
}

static float* rot270(float* in, int & w, int & h, int c) {
  return rot180(rot90(in, w, h, c), w, h, c);
}

//TODO: isolate the dither from the pixel conversion and make adjustable bit depth
static uint8_t* dither8(const float* buffer, int w, int h, int c) {
  uint8_t* ret = (uint8_t*)malloc(w * h * c);
  float* errLine0 = (float*)calloc(w * c, sizeof(float));
  float* errLine1 = (float*)calloc(w * c, sizeof(float));
  float* errLines[2] = { errLine0, errLine1 };
  if (!ret || !errLine0 || !errLine1) {
    std::cerr << "out of memory" << std::endl;
    die(__LINE__);
  }

  for (int y = 0; y < h; ++y) {
    float* errIn = errLines[y % 2];
    float* errOut = errLines[(y + 1) % 2];
    for (int x = 0; x < w; ++x) {
      for (int i = 0; i < c; ++i) {
        const float pxlF = buffer[(w * y + x) * c + i] * 255.0f + errIn[x * c + i];
        auto pxlU8 = llrintf(pxlF);
        if (pxlU8 < 0) pxlU8 = 0;
        if (pxlU8 > 255) pxlU8 = 255;
        ret[(w * y + x) * c + i] = (uint8_t)pxlU8;

        float err = pxlU8 - pxlF;

        if (x < w - 1) errIn[x * c + i] += err * (7.0f / 16.0f);
        if (x > 0) errOut[(x - 1) * c + i] += err * (3.0f / 16.0f);
        errOut[x * c + i] += err * (5.0f / 16.0f);
        if (x < w - 1) errOut[(x + 1) * c + i] += err * (1.0f / 16.0f);
      }
    }
    memset(errIn, 0x00, sizeof(float) * w * c);
  }

  free(errLine0);
  free(errLine1);
  return ret;
}

static uint16_t* dither16(const float* buffer, int w, int h, int c) {
  uint16_t* ret = (uint16_t*)malloc(w * h * c * 2);
  float* errLine0 = (float*)calloc(w * c, sizeof(float));
  float* errLine1 = (float*)calloc(w * c, sizeof(float));
  float* errLines[2] = { errLine0, errLine1 };
  if (!ret || !errLine0 || !errLine1) {
    std::cerr << "out of memory" << std::endl;
    die(__LINE__);
  }

  for (int y = 0; y < h; ++y) {
    float* errIn = errLines[y % 2];
    float* errOut = errLines[(y + 1) % 2];
    for (int x = 0; x < w; ++x) {
      for (int i = 0; i < c; ++i) {
        const float pxlF = buffer[(w * y + x) * c + i] * 65535.0f + errIn[x * c + i];
        auto pxlU16 = llrintf(pxlF);
        if (pxlU16 < 0) pxlU16 = 0;
        if (pxlU16 > 65535) pxlU16 = 65535;
        ret[(w * y + x) * c + i] = (uint16_t)pxlU16;

        float err = pxlU16 - pxlF;

        if (x < w - 1) errIn[x * c + i] += err * (7.0f / 16.0f);
        if (x > 0) errOut[(x - 1) * c + i] += err * (3.0f / 16.0f);
        errOut[x * c + i] += err * (5.0f / 16.0f);
        if (x < w - 1) errOut[(x + 1) * c + i] += err * (1.0f / 16.0f);
      }
    }
    memset(errIn, 0x00, sizeof(float) * w * c);
  }

  free(errLine0);
  free(errLine1);
  return ret;
}

static uint8_t* floor8(const float* buffer, int w, int h, int c) {
  uint8_t* ret = (uint8_t*)malloc(w * h * c);
  if (!ret) {
    std::cerr << "out of memory" << std::endl;
    die(__LINE__);
  }

  for (int i = 0; i < w * h * c; ++i) {
    auto f = floorf(255.0f * buffer[i]);
    if (f < 0) {
      ret[i] = 0;
    }
    else if (f > 255.0f) {
      ret[i] = 255;
    }
    else {
      ret[i] = (uint8_t)f;
    }
  }

  return ret;
}

static uint8_t* round8(const float* buffer, int w, int h, int c) {
  uint8_t* ret = (uint8_t*)malloc(w * h * c);
  if (!ret) {
    std::cerr << "out of memory" << std::endl;
    die(__LINE__);
  }

  for (int i = 0; i < w * h * c; ++i) {
    auto f = roundf(255.0f * buffer[i]);
    if (f < 0) {
      ret[i] = 0;
    }
    else if (f > 255.0f) {
      ret[i] = 255;
    }
    else {
      ret[i] = (uint8_t)f;
    }
  }

  return ret;
}

static float* expand8to32(const uint8_t* input, int w, int h, int c) {
  float* ret = (float*)malloc(w * h * c * sizeof(float));
  if (!ret) {
    std::cerr << "out of memory" << std::endl;
    die(__LINE__);
  }
  constexpr float scalar = (float)(1.0 / 255.0);
  for (int i = 0; i < w * h * c; ++i) {
    ret[i] = input[i] * scalar;
  }
  return ret;
}

static float* expand16to32(const uint16_t* input, int w, int h, int c) {
  float* ret = (float*)malloc(w * h * c * sizeof(float));
  if (!ret) {
    std::cerr << "out of memory" << std::endl;
    die(__LINE__);
  }
  constexpr float scalar = (float)(1.0 / 65535.0);
  for (int i = 0; i < w * h * c; ++i) {
    ret[i] = input[i] * scalar;
  }
  return ret;
}

static void stats(const float* buf, size_t n, int c) {
  float * pixels = (float*)malloc(sizeof(float) * n * c);
  memcpy(pixels, buf, sizeof(float) * n * c);
  std::cout << "Total pixels:" << n << std::endl;
  std::cout << "Unique pixels:";
  qsort(pixels, n, sizeof(float) * c, fltVecCmp);
  unsigned pcount = 1;
  for (size_t i = 1; i < n; ++i) {
    const float* pixel = pixels + i * c;
    pcount += !memcmp(pixel, pixel - c, sizeof(float) * c);
  }
  std::cout << pcount << std::endl;
  std::cout << "Unique/Total Ratio:" << 100 * pcount / n << "%" << std::endl;

  for (int i = 0; i < c; ++i) {
    std::cout << "ch" << i << ":";
    float fmin = INFINITY;
    float fmax = -INFINITY;
    for (size_t j = i; j < n * c; j += c) {
      float pixel = buf[j];
      fmin = std::min(pixel, fmin);
      fmax = std::max(pixel, fmax);
    }
    std::cout << "[" << fmin << "," << fmax << "]/1" << " ";
    std::cout << "[" << (int)(fmin * 255.0f) << ", " << (int)(255.0f * fmax + .999f) << "] / 255";
    std::cout << std::endl;
  }

  free(pixels);
}


int main(int nArgs, const char* args[])
{
  float* buffer = NULL;

  if (nArgs == 1) {
    usage(args[0]);

    std::cerr << "No arguments detected, trying to process test files" << std::endl;

    static const char* dummyArgs[] = { 
      args[0], "-h", "-idct",
      "-i", "input1.jpg", "-s",
      "-o", "output1.fff.data",
      "-o", "output1.png",
      "-o", "output1.jpg",
      "-o16", "output1.161616.png",
      "-rot90",
      "-o", "output1_r90.jpg",
      "-rot180",
      "-o", "output1_r90_180.jpg",
      "-rot270",
      "-o", "output1_r90_180_270.jpg",
      "-dither8",
      "-o", "output1_r90_180_270_dither.jpg",
      "-rot180",
      "-q", "10",
      "-o", "output1_q10.jpg",
      //FFF2JPG
      //"-i", "input2.fff.data",
      //"-o", "output2.fff.data",
      //"-o", "output2.jpg",
    };
    args = dummyArgs;
    nArgs = sizeof(dummyArgs) / sizeof(dummyArgs[0]);
  }

  int w = -1;
  int h = -1;
  int c = -1;
  int q = 99;

  for (int arg = 1; arg < nArgs; ++arg) {
    bool more = arg + 1 < nArgs;

    if (!strcmp(args[arg], "-h") || !strcmp(args[arg], "--help")) {
      usage(args[0]);
    }
    else if (!strcmp(args[arg], "-w") && more) {
      arg++;
      w = strtol(args[arg], NULL, 10);
      std::cout << "width set to " << w << std::endl;
    }
    else if (!strcmp(args[arg], "-q") && more) {
      arg++;
      q = strtol(args[arg], NULL, 10);
      std::cout << "JPEG encoder quality factor set to " << q << std::endl;
    }
    else if (!strcmp(args[arg], "-rot90")) {
      buffer = rot90(buffer, w, h, c);
    }
    else if (!strcmp(args[arg], "-rot180")) {
      buffer = rot180(buffer, w, h, c);
    }
    else if (!strcmp(args[arg], "-rot270")) {
      buffer = rot270(buffer, w, h, c);
    }
    else if (!strcmp(args[arg], "-dither8")) {
      uint8_t* buf8 = dither8(buffer, w, h, c);
      free(buffer);
      buffer = expand8to32(buf8, w, h, c);
      free(buf8);
    }
    else if (!strcmp(args[arg], "-floor8")) {
      uint8_t* buf8 = floor8(buffer, w, h, c);
      free(buffer);
      buffer = expand8to32(buf8, w, h, c);
      free(buf8);
    }
    else if (!strcmp(args[arg], "-round8")) {
      uint8_t* buf8 = round8(buffer, w, h, c);
      free(buffer);
      buffer = expand8to32(buf8, w, h, c);
      free(buf8);
    }
    else if (!strcmp(args[arg], "-o") || !strcmp(args[arg], "-o16")) {
      int bits = strcmp(args[arg], "-o16") ? -1 : 16;

      if (!more) {
        std::cerr << "missing argument to output flag" << std::endl;
        die(__LINE__);
      }

      arg++;
      if (!buffer) {
        std::cerr << "no buffer loaded" << std::endl;
        die(__LINE__);
      }
      if ((strEndsWith(args[arg], ".f.data") && c != 1) || (strEndsWith(args[arg], ".fff.data") && c != 3)) {
        std::cerr << "cannot write " << c << " channel image to " << args[arg] << std::endl;
        die(__LINE__);
      }
      else if (strEndsWith(args[arg], ".fff.data") || strEndsWith(args[arg], ".f.data")) {
        FILE* fp = fopen(args[arg], "wb");
        if (!fp) {
          std::cerr << "unable to open " << args[arg] << std::endl;
          die(__LINE__);
        }
        size_t nFloats = w * h * c;
        if (nFloats != fwrite(buffer, sizeof(float), nFloats, fp)) {
          std::cerr << "fwrite() failure writing " << args[arg] << std::endl;
          die(__LINE__);
        }
        fclose(fp);
        std::cout << "wrote " << args[arg] << std::endl;
      }
      else if (strEndsWith(args[arg], ".jpg")) {
        if (!stbi_write_jpg(args[arg], w, h, c, buffer, q)) {
          std::cerr << "stbi_write_jpg() failure" << std::endl;
          die(__LINE__);
        }
        std::cout << "wrote " << args[arg] << std::endl;
      }
      else if (strEndsWith(args[arg], ".png")) {
        if (bits == 16) {
          std::cerr << "dithering to 16-bit output" << std::endl;
          auto dithered = dither16(buffer, w, h, c);
          std::cerr << "writing 16-bit PNG to "<< args[arg] << std::endl;
          if (!stbi_write_png16(args[arg], w, h, c, dithered, 0)) {
            std::cerr << "stbi_write_png16() failure" << std::endl;
            die(__LINE__);
          }
          free(dithered);
        }
        else {
          std::cerr << "dithering to 8-bit output" << std::endl;
          auto dithered = dither8(buffer, w, h, c);
          std::cerr << "writing 8-bit PNG to "<< args[arg] << std::endl;
          if (!stbi_write_png(args[arg], w, h, c, dithered, 0)) {
            std::cerr << "stbi_write_png() failure" << std::endl;
            die(__LINE__);
          }
          free(dithered);
        }
        std::cout << "wrote " << args[arg] << std::endl;
      }
      else {
        std::cerr << "do not know how to write: " << args[arg] << std::endl;
        die(__LINE__);
      }
    }
    else if (!strcmp(args[arg], "-i")) {
      if (!more) {
        std::cerr << "missing argument to input flag" << std::endl;
        die(__LINE__);
      }

      free(buffer);
      buffer = NULL;

      arg++;
      if (strEndsWith(args[arg], ".data") && w <= 0) {
        std::cerr << "must set width before loading raw image" << std::endl;
        die(__LINE__);
      }
      else if (strEndsWith(args[arg], ".fff.data") || strEndsWith(args[arg], ".f.data")) {
        //load RGB-FFF
        c = strEndsWith(args[arg], ".fff.data") ? 3 : 1;
        h = 0;

        // --- Reading from a file (mode: "r" for read) ---
        FILE* fp = fopen(args[arg], "rb");
        if (fp == NULL) {
          std::cerr << "Error opening file for reading:" << args[arg] << std::endl;
          die(__LINE__);
        }

        size_t r = w * c;
        size_t sz = 0;
        do {
          while (sizeof(float) * w * c * (h + 1) > sz) {
            sz = sz * 2 + sizeof(float) * w * c;
            buffer = (float*)realloc(buffer, sz);
            if (!buffer) {
              std::cerr << "Out of memory" << std::endl;
              die(__LINE__);
            }
          }
          r = fread(buffer + w * c * h, sizeof(float), w * c, fp);
          h++;
        } while (r == (size_t)w * c);
        if (r != 0) {
          std::cerr << "Read " << r << " samples, which is not divisible into " << w << "*" << c << "samples per row" << std::endl;
          std::cerr << "Please verify size(" << args[arg] << ") = width * height * channels * 4B" << std::endl;
          die(__LINE__);
        }
        h--;
        std::cout << "Loaded " << args[arg] << " " << w << "x" << h << " c=" << c << std::endl;
      }
      else if (strEndsWith(args[arg], ".jpg")) {
        std::cout << "Checking channel count:";
        if (!stbi_info(args[arg], &w, &h, &c)) {
          std::cerr << "stbi_info_from_file() failure: " << stbi_failure_reason() << " " << args[arg] << std::endl;
          die(__LINE__);
        }
        std::cout << c << std::endl;

        std::cout << "Reading data from " << args[arg] << std::endl;
        buffer = stbi_loadf(args[arg], &w, &h, &c, c);
        if (!buffer) {
          std::cerr << "stbi_loadf() failure: " << stbi_failure_reason() << " " << args[arg] << std::endl;
          die(__LINE__);
        }
        std::cout << "Loaded " << args[arg] << " " << w << "x"<< h << " c=" << c << std::endl;
      }
      else if (strEndsWith(args[arg], ".png")) {
        std::cout << "Checking channel count:";
        if (!stbi_info(args[arg], &w, &h, &c)) {
          std::cerr << "stbi_info_from_file() failure: " << stbi_failure_reason() << " " << args[arg] << std::endl;
          die(__LINE__);
        }
        std::cout << c << std::endl;

        std::cout << "Reading data from " << args[arg] << std::endl;
        uint16_t* buf16 = stbi_load_16(args[arg], &w, &h, &c, c);
        if (!buf16) {
          std::cerr << "stbi_load_16() failure: " << stbi_failure_reason() << " " << args[arg] << std::endl;
          die(__LINE__);
        }
        buffer = expand16to32(buf16, w, h, c);
        if (!buffer) {
          std::cerr << "expand16to32 failure" << std::endl;
          die(__LINE__);
        }
        free(buf16);
        std::cout << "Loaded " << args[arg] << " " << w << "x" << h << " c=" << c << std::endl;
      }
      else {
        std::cerr << "don't know how to load:" << args[arg] << std::endl;
        die(__LINE__);
      }
    }
    else if (!strcmp(args[arg], "-s")) {
      if (!buffer) {
        std::cerr << "Image buffer empty" << std::endl;
        die(__LINE__);
      }
      stats(buffer, w * h, c);
    }
    else if (!strcmp(args[arg], "-idct")) {
      idct_test();
    }
    else {
      std::cerr << "did not understand arg: " << args[arg] << std::endl;
      die(__LINE__);
    }
  }

  free(buffer);
  buffer = NULL;
  return 0;
}

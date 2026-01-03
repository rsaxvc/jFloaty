#ifndef STBI_INCLUDE_STB_IMAGE_F32_H
#define STBI_INCLUDE_STB_IMAGE_F32_H

#if defined(STBI_NO_JPEG) && defined(STBI_NO_PNG) && defined(STBI_NO_BMP) && defined(STBI_NO_PSD) && defined(STBI_NO_TGA) && defined(STBI_NO_GIF) && defined(STBI_NO_PIC) && defined(STBI_NO_PNM)
// nothing
#else
//////////////////////////////////////////////////////////////////////////////
//
//  generic converter from built-in img_n to req_comp
//    individual types do this automatically as much as possible (e.g. jpeg
//    does all cases internally since it needs to colorspace convert anyway,
//    and it never has alpha, so very few cases ). png can automatically
//    interleave an alpha=255 channel, but falls back to this for other cases
//
//  assume data buffer is malloced, so malloc a new one and free that one
//  only failure mode is malloc failing

static float stbi__compute_y_f32(float r, float g, float b)
{
    return r * (77.0f/256.0f) + g * (150.0f/256.0f) + b * (29.0f/256.0f);
}
#endif
//////////////////////////////////////////////////////////////////////////////
//
//  "baseline" JPEG/JFIF decoder
//
//    simple implementation
//      - doesn't support delayed output of y-dimension
//      - simple interface (only one output format: 8-bit interleaved RGB)
//      - doesn't try to recover corrupt jpegs
//      - doesn't allow partial loading, loading multiple at once
//      - still fast on x86 (copying globals into locals doesn't help x86)
//      - allocates lots of intermediate memory (full size of all components)
//        - non-interleaved case requires this anyway
//        - allows good upsampling (see next)
//    high-quality
//      - upsampled channels are bilinearly interpolated, even across blocks
//      - quality integer IDCT derived from IJG's 'slow'
//    performance
//      - fast huffman; reasonable integer IDCT
//      - some SIMD kernels for common paths on targets with SSE2/NEON
//      - uses a lot of intermediate memory, could cache poorly

#ifndef STBI_NO_JPEG

// derived from jidctint -- DCT_ISLOW
#define STBI__IDCT_1D_F32(s0,s1,s2,s3,s4,s5,s6,s7) \
   float t0,t1,t2,t3,p1,p2,p3,p4,p5,x0,x1,x2,x3; \
   p2 = s2;                                    \
   p3 = s6;                                    \
   p1 = (p2+p3) * (0.5411961f);       \
   t2 = p1 + p3*(-1.847759065f);      \
   t3 = p1 + p2*( 0.765366865f);      \
   p2 = s0;                                    \
   p3 = s4;                                    \
   t0 = (p2+p3);                      \
   t1 = (p2-p3);                      \
   x0 = t0+t3;                                 \
   x3 = t0-t3;                                 \
   x1 = t1+t2;                                 \
   x2 = t1-t2;                                 \
   t0 = s7;                                    \
   t1 = s5;                                    \
   t2 = s3;                                    \
   t3 = s1;                                    \
   p3 = t0+t2;                                 \
   p4 = t1+t3;                                 \
   p1 = t0+t3;                                 \
   p2 = t1+t2;                                 \
   p5 = (p3+p4)*( 1.175875602f);      \
   t0 = t0*( 0.298631336f);           \
   t1 = t1*( 2.053119869f);           \
   t2 = t2*( 3.072711026f);           \
   t3 = t3*( 1.501321110f);           \
   p1 = p5 + p1*(-0.899976223f);      \
   p2 = p5 + p2*(-2.562915447f);      \
   p3 = p3*(-1.961570560f);           \
   p4 = p4*(-0.390180644f);           \
   t3 += p1+p4;                                \
   t2 += p2+p3;                                \
   t1 += p2+p4;                                \
   t0 += p1+p3;

stbi_inline static float stbi__clamp_f32(float x) {
    if (x > 1.0f) return 1.0f;
    if (x < 0.0f) return 0.0f;
    return x;
}

static void stbi__idct_block_f32(float* out, int out_stride, short data[64])
{
  static const float div1020 = 1.0 / 1020.0;
  int i;
  float val[64], * v = val;
  float* o;
  short* d = data;

    // columns
    for (i = 0; i < 8; ++i, ++d, ++v) {
        // if all zeroes, shortcut -- this avoids dequantizing 0s and IDCTing
        if (d[8] == 0 && d[16] == 0 && d[24] == 0 && d[32] == 0
            && d[40] == 0 && d[48] == 0 && d[56] == 0) {
            //    no shortcut                 0     seconds
            //    (1|2|3|4|5|6|7)==0          0     seconds
            //    all separate               -0.047 seconds
            //    1 && 2|3 && 4|5 && 6|7:    -0.047 seconds
            int dcterm = d[0];
            //TODO: test
            v[0] = v[8] = v[16] = v[24] = v[32] = v[40] = v[48] = v[56] = dcterm * div1020;
        } else {
            STBI__IDCT_1D_F32(d[0], d[8], d[16], d[24], d[32], d[40], d[48], d[56])
            // Inputs were fixed point; let's bring them back down
            v[ 0] = (x0 + t3) * div1020;
            v[56] = (x0 - t3) * div1020;
            v[ 8] = (x1 + t2) * div1020;
            v[48] = (x1 - t2) * div1020;
            v[16] = (x2 + t1) * div1020;
            v[40] = (x2 - t1) * div1020;
            v[24] = (x3 + t0) * div1020;
            v[32] = (x3 - t0) * div1020;
        }
    }

    for (i = 0, v = val, o = out; i < 8; ++i, v += 8, o += out_stride) {
        // no fast case since the first 1D IDCT spread components out
        STBI__IDCT_1D_F32(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7])

        o[0] = stbi__clamp_f32((x0 + t3)*.5f+.5f);
        o[7] = stbi__clamp_f32((x0 - t3)*.5f+.5f);
        o[1] = stbi__clamp_f32((x1 + t2)*.5f+.5f);
        o[6] = stbi__clamp_f32((x1 - t2)*.5f+.5f);
        o[2] = stbi__clamp_f32((x2 + t1)*.5f+.5f);
        o[5] = stbi__clamp_f32((x2 - t1)*.5f+.5f);
        o[3] = stbi__clamp_f32((x3 + t0)*.5f+.5f);
        o[4] = stbi__clamp_f32((x3 - t0)*.5f+.5f);
    }                                         
}

static int stbi__parse_entropy_coded_data_f32(stbi__jpeg* z)
{
    stbi__jpeg_reset(z);
    if (!z->progressive) {
        if (z->scan_n == 1) {
            int i, j;
            STBI_SIMD_ALIGN(short, data[64]);
            int n = z->order[0];
            // non-interleaved data, we just need to process one block at a time,
            // in trivial scanline order
            // number of blocks to do just depends on how many actual "pixels" this
            // component has, independent of interleaved MCU blocking and such
            int w = (z->img_comp[n].x + 7) >> 3;
            int h = (z->img_comp[n].y + 7) >> 3;
            for (j = 0; j < h; ++j) {
                for (i = 0; i < w; ++i) {
                    int ha = z->img_comp[n].ha;
                    if (!stbi__jpeg_decode_block(z, data, z->huff_dc + z->img_comp[n].hd, z->huff_ac + ha, z->fast_ac[ha], n, z->dequant[z->img_comp[n].tq])) return 0;
                    z->idct_block_kernel(z->img_comp[n].data + (z->img_comp[n].w2 * j * 8 + i * 8 ) * sizeof(float), z->img_comp[n].w2, data);
                    // every data block is an MCU, so countdown the restart interval
                    if (--z->todo <= 0) {
                        if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                        // if it's NOT a restart, then just bail, so we get corrupt data
                        // rather than no data
                        if (!STBI__RESTART(z->marker)) return 1;
                        stbi__jpeg_reset(z);
                    }
                }
            }
            return 1;
        }
        else { // interleaved
            int i, j, k, x, y;
            STBI_SIMD_ALIGN(short, data[64]);
            for (j = 0; j < z->img_mcu_y; ++j) {
                for (i = 0; i < z->img_mcu_x; ++i) {
                    // scan an interleaved mcu... process scan_n components in order
                    for (k = 0; k < z->scan_n; ++k) {
                        int n = z->order[k];
                        // scan out an mcu's worth of this component; that's just determined
                        // by the basic H and V specified for the component
                        for (y = 0; y < z->img_comp[n].v; ++y) {
                            for (x = 0; x < z->img_comp[n].h; ++x) {
                                int x2 = (i * z->img_comp[n].h + x) * 8;
                                int y2 = (j * z->img_comp[n].v + y) * 8;
                                int ha = z->img_comp[n].ha;
                                if (!stbi__jpeg_decode_block(z, data, z->huff_dc + z->img_comp[n].hd, z->huff_ac + ha, z->fast_ac[ha], n, z->dequant[z->img_comp[n].tq])) return 0;
                                z->idct_block_kernel(z->img_comp[n].data + (z->img_comp[n].w2 * y2 + x2) * sizeof(float), z->img_comp[n].w2, data);
                            }
                        }
                    }
                    // after all interleaved components, that's an interleaved MCU,
                    // so now count down the restart interval
                    if (--z->todo <= 0) {
                        if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                        if (!STBI__RESTART(z->marker)) return 1;
                        stbi__jpeg_reset(z);
                    }
                }
            }
            return 1;
        }
    }
    else {
        if (z->scan_n == 1) {
            int i, j;
            int n = z->order[0];
            // non-interleaved data, we just need to process one block at a time,
            // in trivial scanline order
            // number of blocks to do just depends on how many actual "pixels" this
            // component has, independent of interleaved MCU blocking and such
            int w = (z->img_comp[n].x + 7) >> 3;
            int h = (z->img_comp[n].y + 7) >> 3;
            for (j = 0; j < h; ++j) {
                for (i = 0; i < w; ++i) {
                    short* data = z->img_comp[n].coeff + 64 * (i + j * z->img_comp[n].coeff_w);
                    if (z->spec_start == 0) {
                        if (!stbi__jpeg_decode_block_prog_dc(z, data, &z->huff_dc[z->img_comp[n].hd], n))
                            return 0;
                    }
                    else {
                        int ha = z->img_comp[n].ha;
                        if (!stbi__jpeg_decode_block_prog_ac(z, data, &z->huff_ac[ha], z->fast_ac[ha]))
                            return 0;
                    }
                    // every data block is an MCU, so countdown the restart interval
                    if (--z->todo <= 0) {
                        if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                        if (!STBI__RESTART(z->marker)) return 1;
                        stbi__jpeg_reset(z);
                    }
                }
            }
            return 1;
        }
        else { // interleaved
            int i, j, k, x, y;
            for (j = 0; j < z->img_mcu_y; ++j) {
                for (i = 0; i < z->img_mcu_x; ++i) {
                    // scan an interleaved mcu... process scan_n components in order
                    for (k = 0; k < z->scan_n; ++k) {
                        int n = z->order[k];
                        // scan out an mcu's worth of this component; that's just determined
                        // by the basic H and V specified for the component
                        for (y = 0; y < z->img_comp[n].v; ++y) {
                            for (x = 0; x < z->img_comp[n].h; ++x) {
                                int x2 = (i * z->img_comp[n].h + x);
                                int y2 = (j * z->img_comp[n].v + y);
                                short* data = z->img_comp[n].coeff + 64 * (x2 + y2 * z->img_comp[n].coeff_w);
                                if (!stbi__jpeg_decode_block_prog_dc(z, data, &z->huff_dc[z->img_comp[n].hd], n))
                                    return 0;
                            }
                        }
                    }
                    // after all interleaved components, that's an interleaved MCU,
                    // so now count down the restart interval
                    if (--z->todo <= 0) {
                        if (z->code_bits < 24) stbi__grow_buffer_unsafe(z);
                        if (!STBI__RESTART(z->marker)) return 1;
                        stbi__jpeg_reset(z);
                    }
                }
            }
            return 1;
        }
    }
}

// decode image to YCbCr format
static int stbi__decode_jpeg_image_f32(stbi__jpeg* j)
{
    int m;
    for (m = 0; m < 4; m++) {
        j->img_comp[m].raw_data = NULL;
        j->img_comp[m].raw_coeff = NULL;
    }
    j->restart_interval = 0;
    if (!stbi__decode_jpeg_header(j, STBI__SCAN_load)) return 0;
    m = stbi__get_marker(j);
    while (!stbi__EOI(m)) {
        if (stbi__SOS(m)) {
            if (!stbi__process_scan_header(j)) return 0;
            if (!stbi__parse_entropy_coded_data_f32(j)) return 0;
            if (j->marker == STBI__MARKER_none) {
                j->marker = stbi__skip_jpeg_junk_at_end(j);
                // if we reach eof without hitting a marker, stbi__get_marker() below will fail and we'll eventually return 0
            }
            m = stbi__get_marker(j);
            if (STBI__RESTART(m))
                m = stbi__get_marker(j);
        }
        else if (stbi__DNL(m)) {
            int Ld = stbi__get16be(j->s);
            stbi__uint32 NL = stbi__get16be(j->s);
            if (Ld != 4) return stbi__err("bad DNL len", "Corrupt JPEG");
            if (NL != j->s->img_y) return stbi__err("bad DNL height", "Corrupt JPEG");
            m = stbi__get_marker(j);
        }
        else {
            if (!stbi__process_marker(j, m)) return 1;
            m = stbi__get_marker(j);
        }
    }
    if (j->progressive)
        stbi__jpeg_finish(j);
    return 1;
}


typedef float* (*resample_row_func_f32)(float* out, float* in0, float* in1,
    int w, int hs);

static float* resample_row_1_f32(float* out, float* in_near, float* in_far, int w, int hs)
{
    STBI_NOTUSED(out);
    STBI_NOTUSED(in_far);
    STBI_NOTUSED(w);
    STBI_NOTUSED(hs);
    return in_near;
}

static float* stbi__resample_row_v_2_f32(float* out, float* in_near, float* in_far, int w, int hs)
{
    // need to generate two samples vertically for every one in input
    int i;
    STBI_NOTUSED(hs);
    for (i = 0; i < w; ++i)
        out[i] = .75f * in_near[i] + .25f * in_far[i];
    return out;
}

static float* stbi__resample_row_h_2_f32(float* out, float* in_near, float* in_far, int w, int hs)
{
    // need to generate two samples horizontally for every one in input
    int i;
    float* input = in_near;

    if (w == 1) {
        // if only one sample, can't do any interpolation
        out[0] = out[1] = input[0];
        return out;
    }

    out[0] = input[0];
    out[1] = input[0] * .75f + input[1] * .25f;
    for (i = 1; i < w - 1; ++i) {
        float n = .75f * input[i];
        out[i * 2 + 0] = n + .25f * input[i - 1];
        out[i * 2 + 1] = n + .25f * input[i + 1];
    }
    out[i * 2 + 0] = input[w - 2] * .75f + .25f * input[w - 1];
    out[i * 2 + 1] = input[w - 1];

    STBI_NOTUSED(in_far);
    STBI_NOTUSED(hs);

    return out;
}

static float* stbi__resample_row_hv_2_f32(float* out, float* in_near, float* in_far, int w, int hs)
{  
    // need to generate 2x2 samples for every one in input
    int i;
    float t0, t1;
    if (w == 1) {
        out[0] = out[1] = .75f * in_near[0] + .25f * in_far[0];
        return out;
    }

    out[0] = t1 = .75f * in_near[0] + .25f * in_far[0];
    for (i = 1; i < w; ++i) {
        t0 = t1;
        t1 = .75f * in_near[i] + .25f * in_far[i];
        out[i * 2 - 1] = .75f * t0 + .25f * t1;
        out[i * 2] = .75f * t1 + .25f * t0;
    }
    out[w * 2 - 1] = .25f * t1 + .5f;

    STBI_NOTUSED(hs);

    return out;
}

static float* stbi__resample_row_generic_f32(float* out, float* in_near, float* in_far, int w, int hs)
{
    // resample with nearest-neighbor
    int i, j;
    STBI_NOTUSED(in_far);
    for (i = 0; i < w; ++i)
        for (j = 0; j < hs; ++j)
            out[i * hs + j] = in_near[i];
    return out;
}

static void stbi__YCbCr_to_RGB_row_f32(float* out, const float* y, const float* pcb, const float* pcr, int count, int step)
{
    int i;
    for (i = 0; i < count; ++i) {
        float cr = pcr[i] - .5f;
        float cb = pcb[i] - .5f;
        float r = y[i] + cr * 1.40200f;
        float g = y[i] + cr * -0.71414f + cb * -0.34414f;
        float b = y[i] + cb * 1.77200f;
        r = stbi__clamp_f32(r);
        g = stbi__clamp_f32(g);
        b = stbi__clamp_f32(b);
        out[0] = r;
        out[1] = g;
        out[2] = b;
        out[3] = 1.0f;
        out += step;
    }
}

// set up the kernels
static void stbi__setup_jpeg_f32(stbi__jpeg* j)
{
    j->idct_block_kernel = 
        (void (*)(stbi_uc * out, int out_stride, short data[64]))
        stbi__idct_block_f32;
    j->YCbCr_to_RGB_kernel = 
        (void (*)(stbi_uc * out, const stbi_uc * y, const stbi_uc * pcb, const stbi_uc * pcr, int count, int step))
        stbi__YCbCr_to_RGB_row_f32;
    j->resample_row_hv_2_kernel =
        (stbi_uc * (*)(stbi_uc * out, stbi_uc * in_near, stbi_uc * in_far, int w, int hs))
        stbi__resample_row_hv_2_f32;
}

typedef struct
{
    resample_row_func_f32 resample;
    float* line0, * line1;
    int hs, vs;   // expansion factor in each axis
    int w_lores; // horizontal pixels pre-expansion
    int ystep;   // how far through vertical expansion we are
    int ypos;    // which pre-expansion row we're on
} stbi__resample_f32;

static void f32set(float* d, float v, size_t n) {
  while (n) {
    *d++ = v;
    n--;
  }
}

static float* load_jpeg_image_f32(stbi__jpeg* z, int* out_x, int* out_y, int* comp, int req_comp)
{
    int n, decode_n, is_rgb;
    z->s->img_n = 0; // make stbi__cleanup_jpeg safe

    // validate req_comp
    if (req_comp < 0 || req_comp > 4) return (float*)stbi__errpuc("bad req_comp", "Internal error");

    // load a jpeg image from whichever source, but leave in YCbCr format
    if (!stbi__decode_jpeg_image_f32(z)) { stbi__cleanup_jpeg(z); return NULL; }

    // determine actual number of components to generate
    n = req_comp ? req_comp : z->s->img_n >= 3 ? 3 : 1;

    is_rgb = z->s->img_n == 3 && (z->rgb == 3 || (z->app14_color_transform == 0 && !z->jfif));

    if (z->s->img_n == 3 && n < 3 && !is_rgb)
        decode_n = 1;
    else
        decode_n = z->s->img_n;

    // nothing to do if no components requested; check this now to avoid
    // accessing uninitialized coutput[0] later
    if (decode_n <= 0) { stbi__cleanup_jpeg(z); return NULL; }

    // resample and color-convert
    {
        int k;
        unsigned int i, j;
        float* output;
        float* coutput[4] = { NULL, NULL, NULL, NULL };

        stbi__resample_f32 res_comp[4];

        for (k = 0; k < decode_n; ++k) {
            stbi__resample_f32* r = &res_comp[k];

            // allocate line buffer big enough for upsampling off the edges
            // with upsample factor of 4
            z->img_comp[k].linebuf = (stbi_uc*)stbi__malloc_mad2(sizeof(float), z->s->img_x, 3);
            if (!z->img_comp[k].linebuf) { stbi__cleanup_jpeg(z); return (float*)stbi__errpuc("outofmem", "Out of memory"); }

            r->hs = z->img_h_max / z->img_comp[k].h;
            r->vs = z->img_v_max / z->img_comp[k].v;
            r->ystep = r->vs >> 1;
            r->w_lores = (z->s->img_x + r->hs - 1) / r->hs;
            r->ypos = 0;
            r->line0 = r->line1 = (float*)(z->img_comp[k].data);

            if (r->hs == 1 && r->vs == 1) r->resample = resample_row_1_f32;
            else if (r->hs == 1 && r->vs == 2) r->resample = stbi__resample_row_v_2_f32;
            else if (r->hs == 2 && r->vs == 1) r->resample = stbi__resample_row_h_2_f32;
            else if (r->hs == 2 && r->vs == 2) r->resample = (resample_row_func_f32)z->resample_row_hv_2_kernel;
            else                               r->resample = stbi__resample_row_generic_f32;
        }

        // can't error after this so, this is safe
        output = (float*)stbi__malloc_mad4(sizeof(float), n, z->s->img_x, z->s->img_y, 1);
        if (!output) { stbi__cleanup_jpeg(z); return (float*)stbi__errpuc("outofmem", "Out of memory"); }

        // now go ahead and resample
        for (j = 0; j < z->s->img_y; ++j) {
            float* out = output + n * z->s->img_x * j;
            for (k = 0; k < decode_n; ++k) {
                stbi__resample_f32* r = &res_comp[k];
                int y_bot = r->ystep >= (r->vs >> 1);
                coutput[k] = r->resample((float*)(z->img_comp[k].linebuf),
                    y_bot ? r->line1 : r->line0,
                    y_bot ? r->line0 : r->line1,
                    r->w_lores, r->hs);
                if (++r->ystep >= r->vs) {
                    r->ystep = 0;
                    r->line0 = r->line1;
                    if (++r->ypos < z->img_comp[k].y)
                        r->line1 += z->img_comp[k].w2;
                }
            }
            if (n >= 3) {
                float * y = coutput[0];
                if (z->s->img_n == 3) {
                    if (is_rgb) {
                        for (i = 0; i < z->s->img_x; ++i) {
                            out[0] = y[i];
                            out[1] = coutput[1][i];
                            out[2] = coutput[2][i];
                            out[3] = 1.0;
                            out += n;
                        }
                    }
                    else {
                      f32set(out, j % 2, z->s->img_x * n);
                      z->YCbCr_to_RGB_kernel((stbi_uc*)out, (stbi_uc*)y, (stbi_uc*)coutput[1], (stbi_uc*)coutput[2], z->s->img_x, n);
                    }
                }
                else if (z->s->img_n == 4) {
                    if (z->app14_color_transform == 0) { // CMYK
                        for (i = 0; i < z->s->img_x; ++i) {
                            float m = coutput[3][i];
                            out[0] = coutput[0][i] * m;
                            out[1] = coutput[1][i] * m;
                            out[2] = coutput[2][i] * m;
                            out[3] = 1.0f;
                            out += n;
                        }
                    }
                    else if (z->app14_color_transform == 2) { // YCCK
                        z->YCbCr_to_RGB_kernel((stbi_uc*)out, (stbi_uc*)y, (stbi_uc*)coutput[1], (stbi_uc*)coutput[2], z->s->img_x, n);
                        for (i = 0; i < z->s->img_x; ++i) {
                            float m = coutput[3][i];
                            out[0] = (1.0f - out[0]) * m;
                            out[1] = (1.0f - out[1]) * m;
                            out[2] = (1.0f - out[2]) * m;
                            out += n;
                        }
                    }
                    else { // YCbCr + alpha?  Ignore the fourth channel for now
                        z->YCbCr_to_RGB_kernel((stbi_uc*)out, (stbi_uc*)y, (stbi_uc*)coutput[1], (stbi_uc*)coutput[2], z->s->img_x, n);
                    }
                }
                else
                    for (i = 0; i < z->s->img_x; ++i) {
                        out[0] = out[1] = out[2] = y[i];
                        out[3] = 1.0f; // not used if n==3
                        out += n;
                    }
            }
            else {
                if (is_rgb) {
                    if (n == 1)
                        for (i = 0; i < z->s->img_x; ++i)
                            *out++ = stbi__compute_y_f32(coutput[0][i], coutput[1][i], coutput[2][i]);
                    else {
                        for (i = 0; i < z->s->img_x; ++i, out += 2) {
                            out[0] = stbi__compute_y_f32(coutput[0][i], coutput[1][i], coutput[2][i]);
                            out[1] = 1.0f;
                        }
                    }
                }
                else if (z->s->img_n == 4 && z->app14_color_transform == 0) {
                    for (i = 0; i < z->s->img_x; ++i) {
                        float m = coutput[3][i];
                        float r = coutput[0][i] * m;
                        float g = coutput[1][i] * m;
                        float b = coutput[2][i] * m;
                        out[0] = stbi__compute_y_f32(r, g, b);
                        out[1] = 1.0f;
                        out += n;
                    }
                }
                else if (z->s->img_n == 4 && z->app14_color_transform == 2) {
                    for (i = 0; i < z->s->img_x; ++i) {
                        out[0] = (1.0f - coutput[0][i]) * coutput[3][i];
                        out[1] = 1.0f;
                        out += n;
                    }
                }
                else {
                    float* y = coutput[0];
                    if (n == 1)
                        for (i = 0; i < z->s->img_x; ++i) 
                            out[i] = y[i];
                    else
                        for (i = 0; i < z->s->img_x; ++i)
                            { *out++ = y[i]; *out++ = 1.0f; }
                }
            }
        }
        stbi__cleanup_jpeg(z);
        *out_x = z->s->img_x;
        *out_y = z->s->img_y;
        if (comp) *comp = z->s->img_n >= 3 ? 3 : 1; // report original components, not output
        return output;
    }
}

static void* stbi__jpeg_load_f32(stbi__context* s, int* x, int* y, int* comp, int req_comp, stbi__result_info* ri)
{
    stbi__jpeg* j = (stbi__jpeg*)stbi__malloc(sizeof(stbi__jpeg));
    if (!j) return stbi__errpuc("outofmem", "Out of memory");
    memset(j, 0, sizeof(stbi__jpeg));
    STBI_NOTUSED(ri);
    j->s = s;
    stbi__setup_jpeg_f32(j);
    void * result = load_jpeg_image_f32(j, x, y, comp, req_comp);
    STBI_FREE(j);
    return result;
}

#endif

#endif //STBI_INCLUDE_STB_IMAGE_F32_H
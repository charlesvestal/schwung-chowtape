/*
 * CHOWTape - Analog tape model audio effect for Move Anything
 *
 * Pure C port of the AnalogTapeModel by Jatin Chowdhury
 * https://github.com/jatinchowdhury18/AnalogTapeModel
 *
 * Core DSP: Jiles-Atherton hysteresis model with RK4 solver,
 * tape loss filter, wow & flutter, chew, and degradation.
 *
 * Original code: Copyright (c) 2019 Jatin Chowdhury (GPL-3.0)
 * This port: Copyright (c) 2026 Charles Vestal (GPL-3.0)
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <stdint.h>
#include "plugin_api_v1.h"

/* ---- Constants ---- */

#define SAMPLE_RATE 44100.0
#define kTwoPi 6.283185307179586
#define kPi 3.141592653589793
#define ONE_THIRD (1.0 / 3.0)
#define NEG_TWO_OVER_15 (-2.0 / 15.0)

#define AUDIO_FX_API_VERSION_2 2
#define AUDIO_FX_INIT_V2_SYMBOL "move_audio_fx_init_v2"

/* Wow/flutter delay line */
#define WF_MAX_DELAY_SAMPLES (1 << 17)  /* ~3 seconds at 44100 */
#define WF_DELAY_MASK (WF_MAX_DELAY_SAMPLES - 1)

/* ---- Host API ---- */

static const host_api_v1_t *g_host = NULL;

static void host_log(const char *msg) {
    if (g_host && g_host->log) g_host->log(msg);
}

/* ---- Simple RNG (LCG) ---- */

static inline uint32_t lcg_next(uint32_t *state) {
    *state = *state * 1664525u + 1013904223u;
    return *state;
}

static inline float lcg_float(uint32_t *state) {
    return (float)(lcg_next(state) & 0x7FFFFF) / (float)0x7FFFFF;
}

/* ---- Smoothed Value ---- */

typedef struct {
    float current;
    float target;
    float step;
    int steps_remaining;
    int total_steps;
} smooth_val_t;

static inline void smooth_init(smooth_val_t *s, float val, int steps) {
    s->current = val;
    s->target = val;
    s->step = 0.0f;
    s->steps_remaining = 0;
    s->total_steps = steps;
}

static inline void smooth_set(smooth_val_t *s, float target) {
    if (target == s->target) return;
    s->target = target;
    s->steps_remaining = s->total_steps;
    s->step = (target - s->current) / (float)s->total_steps;
}

static inline float smooth_next(smooth_val_t *s) {
    if (s->steps_remaining > 0) {
        s->current += s->step;
        s->steps_remaining--;
        if (s->steps_remaining == 0) s->current = s->target;
    }
    return s->current;
}

/* ==================================================================
 * Hysteresis - Jiles-Atherton model with RK4 solver
 * Based on: https://ccrma.stanford.edu/~jatin/420/tape/TapeModel_DAFx.pdf
 * ================================================================== */

typedef struct {
    /* Jiles-Atherton parameters */
    double M_s;       /* Saturation magnetization */
    double a;         /* Shape parameter */
    double k;         /* Pinning parameter */
    double c;         /* Reversibility */
    double nc;        /* 1 - c */

    /* Pre-computed derived values */
    double M_s_oa;
    double M_s_oa_talpha;
    double M_s_oa_tc;
    double M_s_oa_tc_talpha;
    double M_s_oaSq_tc_talpha;
    double M_s_oaSq_tc_talphaSq;

    /* State for hysteresisFunc (cached between calls) */
    double Q, M_diff, L_prime, kap1, f1Denom, f1, f2, f3;
    double coth;
    int nearZero;

    /* Solver state */
    double M_n1;       /* Previous magnetization */
    double H_n1;       /* Previous magnetic field */
    double H_d_n1;     /* Previous field derivative */

    /* Timing */
    double T;          /* 1/sampleRate */
    double Talpha;     /* T/1.9 (for NR solver) */
    double upperLim;
} hysteresis_t;

static const double HYST_ALPHA = 1.6e-3;

static inline int hyst_sign(double x) {
    return (x > 0.0) - (x < 0.0);
}

static inline double hyst_deriv(double x_n, double x_n1, double x_d_n1, double T) {
    const double dAlpha = 0.75;
    return ((1.0 + dAlpha) / T) * (x_n - x_n1) - dAlpha * x_d_n1;
}

static inline double hyst_langevin(hysteresis_t *hp) {
    return !hp->nearZero ? (hp->coth - (1.0 / hp->Q)) : hp->Q * ONE_THIRD;
}

static inline double hyst_langevinD(hysteresis_t *hp) {
    double oneOverQSq = 1.0 / (hp->Q * hp->Q);
    double cothSq = hp->coth * hp->coth;
    return !hp->nearZero ? (oneOverQSq - cothSq + 1.0) : ONE_THIRD;
}

static double hyst_func(double M, double H, double H_d, hysteresis_t *hp) {
    hp->Q = (H + M * HYST_ALPHA) / hp->a;

    if (hp->Q < 0.001 && hp->Q > -0.001) {
        hp->nearZero = 1;
        hp->coth = 0.0;
    } else {
        hp->nearZero = 0;
        hp->coth = 1.0 / tanh(hp->Q);
    }

    double L = hyst_langevin(hp);
    hp->M_diff = L * hp->M_s - M;

    double delta = (H_d >= 0.0) ? 1.0 : -1.0;
    double delta_M = (double)(hyst_sign(delta) == hyst_sign(hp->M_diff));
    hp->kap1 = hp->nc * delta_M;

    hp->L_prime = hyst_langevinD(hp);

    hp->f1Denom = (hp->nc * delta) * hp->k - HYST_ALPHA * hp->M_diff;
    hp->f1 = hp->kap1 * hp->M_diff / hp->f1Denom;
    hp->f2 = hp->L_prime * hp->M_s_oa_tc;
    hp->f3 = 1.0 - (hp->L_prime * hp->M_s_oa_tc_talpha);

    return H_d * (hp->f1 + hp->f2) / hp->f3;
}

static double hyst_rk4(hysteresis_t *hp, double H, double H_d) {
    double T = hp->T;
    double H_1_2 = (H + hp->H_n1) * 0.5;
    double H_d_1_2 = (H_d + hp->H_d_n1) * 0.5;

    double k1 = hyst_func(hp->M_n1, hp->H_n1, hp->H_d_n1, hp) * T;
    double k2 = hyst_func(hp->M_n1 + k1 * 0.5, H_1_2, H_d_1_2, hp) * T;
    double k3 = hyst_func(hp->M_n1 + k2 * 0.5, H_1_2, H_d_1_2, hp) * T;
    double k4 = hyst_func(hp->M_n1 + k3, H, H_d, hp) * T;

    return hp->M_n1 + k1 / 6.0 + k2 / 3.0 + k3 / 3.0 + k4 / 6.0;
}

static inline double hyst_process(hysteresis_t *hp, double H) {
    double H_d = hyst_deriv(H, hp->H_n1, hp->H_d_n1, hp->T);

    double M = hyst_rk4(hp, H, H_d);

    /* Check for instability */
    if (isnan(M) || M > hp->upperLim) {
        M = 0.0;
        H_d = 0.0;
    }

    hp->M_n1 = M;
    hp->H_n1 = H;
    hp->H_d_n1 = H_d;

    return M;
}

static void hyst_cook(hysteresis_t *hp, double drive, double width, double sat) {
    hp->M_s = 0.5 + 1.5 * (1.0 - sat);
    hp->a = hp->M_s / (0.01 + 6.0 * drive);
    hp->c = sqrt(1.0 - width) - 0.01;
    if (hp->c < 0.001) hp->c = 0.001;
    hp->k = 0.47875;
    hp->upperLim = 20.0;

    hp->nc = 1.0 - hp->c;
    hp->M_s_oa = hp->M_s / hp->a;
    hp->M_s_oa_talpha = HYST_ALPHA * hp->M_s_oa;
    hp->M_s_oa_tc = hp->c * hp->M_s_oa;
    hp->M_s_oa_tc_talpha = HYST_ALPHA * hp->M_s_oa_tc;
    hp->M_s_oaSq_tc_talpha = hp->M_s_oa_tc_talpha / hp->a;
    hp->M_s_oaSq_tc_talphaSq = HYST_ALPHA * hp->M_s_oaSq_tc_talpha;
}

static void hyst_init(hysteresis_t *hp, double sr) {
    hp->T = 1.0 / sr;
    hp->Talpha = hp->T / 1.9;
    hp->M_n1 = 0.0;
    hp->H_n1 = 0.0;
    hp->H_d_n1 = 0.0;
    hp->coth = 0.0;
    hp->nearZero = 0;
    hyst_cook(hp, 0.5, 0.5, 0.5);
}

/* ==================================================================
 * Tone Control - Shelf filter (simplified from ToneControl.cpp)
 * ================================================================== */

typedef struct {
    float b0, b1, b2, a1, a2;
    float x1, x2, y1, y2;
} biquad_t;

static void biquad_reset(biquad_t *bq) {
    bq->x1 = bq->x2 = bq->y1 = bq->y2 = 0.0f;
}

static inline float biquad_process(biquad_t *bq, float x) {
    float y = bq->b0 * x + bq->b1 * bq->x1 + bq->b2 * bq->x2
              - bq->a1 * bq->y1 - bq->a2 * bq->y2;
    bq->x2 = bq->x1;
    bq->x1 = x;
    bq->y2 = bq->y1;
    bq->y1 = y;
    return y;
}

/* Low shelf filter */
static void biquad_lowshelf(biquad_t *bq, float fs, float fc, float gain_db) {
    float A = powf(10.0f, gain_db / 40.0f);
    float w0 = (float)(kTwoPi) * fc / fs;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);
    float alpha = sinw0 / 2.0f * sqrtf(2.0f);
    float sqrtA = sqrtf(A);

    float a0 = (A + 1.0f) + (A - 1.0f) * cosw0 + 2.0f * sqrtA * alpha;
    bq->b0 = (A * ((A + 1.0f) - (A - 1.0f) * cosw0 + 2.0f * sqrtA * alpha)) / a0;
    bq->b1 = (2.0f * A * ((A - 1.0f) - (A + 1.0f) * cosw0)) / a0;
    bq->b2 = (A * ((A + 1.0f) - (A - 1.0f) * cosw0 - 2.0f * sqrtA * alpha)) / a0;
    bq->a1 = (-2.0f * ((A - 1.0f) + (A + 1.0f) * cosw0)) / a0;
    bq->a2 = ((A + 1.0f) + (A - 1.0f) * cosw0 - 2.0f * sqrtA * alpha) / a0;
}

/* High shelf filter */
static void biquad_highshelf(biquad_t *bq, float fs, float fc, float gain_db) {
    float A = powf(10.0f, gain_db / 40.0f);
    float w0 = (float)(kTwoPi) * fc / fs;
    float cosw0 = cosf(w0);
    float sinw0 = sinf(w0);
    float alpha = sinw0 / 2.0f * sqrtf(2.0f);
    float sqrtA = sqrtf(A);

    float a0 = (A + 1.0f) - (A - 1.0f) * cosw0 + 2.0f * sqrtA * alpha;
    bq->b0 = (A * ((A + 1.0f) + (A - 1.0f) * cosw0 + 2.0f * sqrtA * alpha)) / a0;
    bq->b1 = (-2.0f * A * ((A - 1.0f) + (A + 1.0f) * cosw0)) / a0;
    bq->b2 = (A * ((A + 1.0f) + (A - 1.0f) * cosw0 - 2.0f * sqrtA * alpha)) / a0;
    bq->a1 = (2.0f * ((A - 1.0f) - (A + 1.0f) * cosw0)) / a0;
    bq->a2 = ((A + 1.0f) - (A - 1.0f) * cosw0 - 2.0f * sqrtA * alpha) / a0;
}

/* ==================================================================
 * Loss Filter - Models head spacing/thickness/gap frequency loss
 * Simplified as a variable-frequency 1-pole lowpass
 * ================================================================== */

typedef struct {
    float z1;
    float a1;
    float b0;
} onepole_t;

static void onepole_reset(onepole_t *f) {
    f->z1 = 0.0f;
    f->a1 = 0.0f;
    f->b0 = 1.0f;
}

static void onepole_set_freq(onepole_t *f, float fc, float fs) {
    float wc = (float)kTwoPi * fc / fs;
    float c = 1.0f / tanf(wc / 2.0f);
    float a0 = c + 1.0f;
    f->b0 = 1.0f / a0;
    f->a1 = (1.0f - c) / a0;
}

static inline float onepole_process(onepole_t *f, float x) {
    float y = f->z1 + x * f->b0;
    f->z1 = x * f->b0 - y * f->a1;
    return y;
}

/* ==================================================================
 * DC Blocker - Simple 1-pole highpass
 * ================================================================== */

typedef struct {
    float x1;
    float y1;
    float R;
} dc_blocker_t;

static void dcblock_init(dc_blocker_t *dc, float fc, float fs) {
    dc->x1 = 0.0f;
    dc->y1 = 0.0f;
    dc->R = 1.0f - ((float)kTwoPi * fc / fs);
    if (dc->R < 0.9f) dc->R = 0.9f;
    if (dc->R > 0.9999f) dc->R = 0.9999f;
}

static inline float dcblock_process(dc_blocker_t *dc, float x) {
    float y = x - dc->x1 + dc->R * dc->y1;
    dc->x1 = x;
    dc->y1 = y;
    return y;
}

/* ==================================================================
 * Wow & Flutter - LFO-modulated delay with Ornstein-Uhlenbeck process
 * ================================================================== */

typedef struct {
    /* Delay buffer (per channel) */
    float buf_l[WF_MAX_DELAY_SAMPLES];
    float buf_r[WF_MAX_DELAY_SAMPLES];
    int write_pos;

    /* Wow LFO state (per channel) */
    float wow_phase[2];
    float wow_angle_delta;
    float wow_amp;

    /* Flutter LFO state */
    float flutter_phase1[2];
    float flutter_phase2[2];
    float flutter_phase3[2];
    float flutter_angle1;
    float flutter_angle2;
    float flutter_angle3;
    float flutter_amp1;
    float flutter_amp2;
    float flutter_amp3;
    float flutter_dc_offset;

    /* Ornstein-Uhlenbeck process for wow variance */
    float oh_y[2];
    float oh_sqrtdelta;
    float oh_T;
    float oh_amt;
    float oh_damping;
    float oh_mean;

    /* Smoothed depths */
    smooth_val_t wow_depth;
    smooth_val_t flutter_depth;

    /* DC blockers */
    dc_blocker_t dc_l;
    dc_blocker_t dc_r;

    float fs;
    uint32_t rng;
} wowflutter_t;

static void wf_init(wowflutter_t *wf, float fs) {
    memset(wf->buf_l, 0, sizeof(wf->buf_l));
    memset(wf->buf_r, 0, sizeof(wf->buf_r));
    wf->write_pos = 0;
    wf->fs = fs;

    wf->wow_phase[0] = wf->wow_phase[1] = 0.0f;
    wf->flutter_phase1[0] = wf->flutter_phase1[1] = 0.0f;
    wf->flutter_phase2[0] = wf->flutter_phase2[1] = 0.0f;
    wf->flutter_phase3[0] = wf->flutter_phase3[1] = 0.0f;

    /* Flutter amplitudes scaled to sample rate */
    wf->flutter_amp1 = -230.0f * 1000.0f / fs;
    wf->flutter_amp2 = -80.0f * 1000.0f / fs;
    wf->flutter_amp3 = -99.0f * 1000.0f / fs;
    wf->flutter_dc_offset = 350.0f * 1000.0f / fs;

    /* Wow amplitude */
    wf->wow_amp = 1000.0f * 1000.0f / fs;

    /* O-U process */
    wf->oh_y[0] = wf->oh_y[1] = 1.0f;
    wf->oh_sqrtdelta = 1.0f / sqrtf(fs);
    wf->oh_T = 1.0f / fs;
    wf->oh_amt = 0.0f;
    wf->oh_damping = 1.0f;
    wf->oh_mean = 0.0f;

    smooth_init(&wf->wow_depth, 0.001f, (int)(fs * 0.05f));
    smooth_init(&wf->flutter_depth, 0.001f, (int)(fs * 0.05f));

    dcblock_init(&wf->dc_l, 15.0f, fs);
    dcblock_init(&wf->dc_r, 15.0f, fs);

    wf->rng = 42;
}

static void wf_set_params(wowflutter_t *wf, float flutter_param) {
    /* Flutter param 0-1 maps to combined wow + flutter */
    float wow_depth_val = powf(flutter_param * 0.7f, 3.0f);
    float flutter_depth_val = powf(powf(flutter_param * 0.8f, 3.0f) * 81.0f / 625.0f, 0.5f);

    smooth_set(&wf->wow_depth, fmaxf(0.001f, wow_depth_val));
    smooth_set(&wf->flutter_depth, fmaxf(0.001f, flutter_depth_val));

    /* Wow rate: ~0.5 to ~4.5 Hz */
    float wow_freq = powf(4.5f, 0.25f) - 1.0f;
    /* Add drift */
    float drift = lcg_float(&wf->rng);
    float freq_adjust = wow_freq * (1.0f + powf(drift, 1.25f) * 0.3f);
    wf->wow_angle_delta = (float)kTwoPi * freq_adjust / wf->fs;

    /* Flutter rate: ~10 Hz */
    float flutter_freq = 0.1f * powf(1000.0f, 0.3f);
    wf->flutter_angle1 = (float)kTwoPi * flutter_freq / wf->fs;
    wf->flutter_angle2 = 2.0f * wf->flutter_angle1;
    wf->flutter_angle3 = 3.0f * wf->flutter_angle1;

    /* O-U process params for wow variance */
    float var = flutter_param * 0.3f;
    float amt_pow = powf(var, 1.25f);
    wf->oh_amt = amt_pow;
    wf->oh_damping = amt_pow * 20.0f + 1.0f;
    wf->oh_mean = amt_pow;
}

/* Cubic Hermite interpolation */
static inline float cubic_interp(float y0, float y1, float y2, float y3, float frac) {
    float c0 = y1;
    float c1 = 0.5f * (y2 - y0);
    float c2 = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
    float c3 = 0.5f * (y3 - y0) + 1.5f * (y1 - y2);
    return ((c3 * frac + c2) * frac + c1) * frac + c0;
}

static float wf_read_delay(float *buf, int write_pos, float delay_samples) {
    float int_part;
    float frac = modff(delay_samples, &int_part);
    int d = (int)int_part;

    int i0 = (write_pos - d - 1) & WF_DELAY_MASK;
    int i1 = (write_pos - d) & WF_DELAY_MASK;
    int i2 = (write_pos - d + 1) & WF_DELAY_MASK;
    int i3 = (write_pos - d + 2) & WF_DELAY_MASK;

    return cubic_interp(buf[i0], buf[i1], buf[i2], buf[i3], frac);
}

static void wf_process_sample(wowflutter_t *wf, float *l, float *r) {
    /* O-U process for wow variance */
    float noise = lcg_float(&wf->rng) * 2.0f - 1.0f;
    float oh_noise = noise / 2.33f;

    for (int ch = 0; ch < 2; ch++) {
        wf->oh_y[ch] += wf->oh_sqrtdelta * oh_noise * wf->oh_amt;
        wf->oh_y[ch] += wf->oh_damping * (wf->oh_mean - wf->oh_y[ch]) * wf->oh_T;
    }

    /* Wow LFO */
    float wow_depth = smooth_next(&wf->wow_depth) * wf->wow_amp;
    float wow_l = wow_depth * (cosf(wf->wow_phase[0]) + wf->oh_y[0] * 0.5f);
    float wow_r = wow_depth * (cosf(wf->wow_phase[1]) + wf->oh_y[1] * 0.5f);

    wf->wow_phase[0] += wf->wow_angle_delta;
    wf->wow_phase[1] += wf->wow_angle_delta;
    if (wf->wow_phase[0] >= (float)kTwoPi) wf->wow_phase[0] -= (float)kTwoPi;
    if (wf->wow_phase[1] >= (float)kTwoPi) wf->wow_phase[1] -= (float)kTwoPi;

    /* Flutter LFO (3 cosine oscillators) */
    static const float phaseOff2 = 13.0f * (float)kPi / 4.0f;
    static const float phaseOff3 = -(float)kPi / 10.0f;

    float flutter_depth = smooth_next(&wf->flutter_depth);
    float flutter_l = flutter_depth * (
        wf->flutter_amp1 * cosf(wf->flutter_phase1[0]) +
        wf->flutter_amp2 * cosf(wf->flutter_phase2[0] + phaseOff2) +
        wf->flutter_amp3 * cosf(wf->flutter_phase3[0] + phaseOff3)
    );
    float flutter_r = flutter_depth * (
        wf->flutter_amp1 * cosf(wf->flutter_phase1[1]) +
        wf->flutter_amp2 * cosf(wf->flutter_phase2[1] + phaseOff2) +
        wf->flutter_amp3 * cosf(wf->flutter_phase3[1] + phaseOff3)
    );

    for (int ch = 0; ch < 2; ch++) {
        wf->flutter_phase1[ch] += wf->flutter_angle1;
        wf->flutter_phase2[ch] += wf->flutter_angle2;
        wf->flutter_phase3[ch] += wf->flutter_angle3;
        if (wf->flutter_phase1[ch] >= (float)kTwoPi) wf->flutter_phase1[ch] -= (float)kTwoPi;
        if (wf->flutter_phase2[ch] >= (float)kTwoPi) wf->flutter_phase2[ch] -= (float)kTwoPi;
        if (wf->flutter_phase3[ch] >= (float)kTwoPi) wf->flutter_phase3[ch] -= (float)kTwoPi;
    }

    /* Total delay in samples */
    float delay_l = (wow_l + flutter_l + flutter_depth * wf->flutter_dc_offset + wow_depth) * wf->fs / 1000.0f;
    float delay_r = (wow_r + flutter_r + flutter_depth * wf->flutter_dc_offset + wow_depth) * wf->fs / 1000.0f;

    if (delay_l < 0.0f) delay_l = 0.0f;
    if (delay_l > (float)(WF_MAX_DELAY_SAMPLES - 4)) delay_l = (float)(WF_MAX_DELAY_SAMPLES - 4);
    if (delay_r < 0.0f) delay_r = 0.0f;
    if (delay_r > (float)(WF_MAX_DELAY_SAMPLES - 4)) delay_r = (float)(WF_MAX_DELAY_SAMPLES - 4);

    /* Write to delay buffers */
    wf->buf_l[wf->write_pos] = *l;
    wf->buf_r[wf->write_pos] = *r;

    /* Read from delay buffers */
    *l = wf_read_delay(wf->buf_l, wf->write_pos, delay_l);
    *r = wf_read_delay(wf->buf_r, wf->write_pos, delay_r);

    wf->write_pos = (wf->write_pos + 1) & WF_DELAY_MASK;

    /* DC blocking */
    *l = dcblock_process(&wf->dc_l, *l);
    *r = dcblock_process(&wf->dc_r, *r);
}

/* ==================================================================
 * Chew - Tape damage simulation (dropout + lowpass)
 * ================================================================== */

typedef struct {
    /* Dropout */
    smooth_val_t dropout_mix;
    smooth_val_t dropout_power;

    /* Lowpass per channel */
    onepole_t filt[2];

    /* Timing */
    int samples_until_change;
    int is_crinkled;
    int sample_counter;
    uint32_t rng;
    float fs;

    /* Params */
    float depth;
    float freq;
    float variance;
} chew_t;

static void chew_init(chew_t *ch, float fs) {
    smooth_init(&ch->dropout_mix, 0.0f, (int)(fs * 0.01f));
    smooth_init(&ch->dropout_power, 1.0f, (int)(fs * 0.005f));
    onepole_reset(&ch->filt[0]);
    onepole_reset(&ch->filt[1]);
    onepole_set_freq(&ch->filt[0], 20000.0f, fs);
    onepole_set_freq(&ch->filt[1], 20000.0f, fs);
    ch->samples_until_change = (int)(fs * 0.5f);
    ch->is_crinkled = 0;
    ch->sample_counter = 0;
    ch->rng = 12345;
    ch->fs = fs;
    ch->depth = 0.0f;
    ch->freq = 0.0f;
    ch->variance = 0.0f;
}

static int chew_get_dry_time(chew_t *ch) {
    float tScale = powf(ch->freq, 0.1f);
    float varScale = powf(lcg_float(&ch->rng) * 2.0f, ch->variance);
    int min_t = (int)((1.0f - tScale) * ch->fs * varScale);
    int max_t = (int)((2.0f - 1.99f * tScale) * ch->fs * varScale);
    if (max_t <= min_t) max_t = min_t + 1;
    return min_t + (int)(lcg_float(&ch->rng) * (float)(max_t - min_t));
}

static int chew_get_wet_time(chew_t *ch) {
    float tScale = powf(ch->freq, 0.1f);
    float start = 0.2f + 0.8f * ch->depth;
    float end = start - (0.001f + 0.01f * ch->depth);
    float varScale = powf(lcg_float(&ch->rng) * 2.0f, ch->variance);
    int min_t = (int)((1.0f - tScale) * ch->fs * varScale);
    int max_t = (int)(((1.0f - tScale) + start - end * tScale) * ch->fs * varScale);
    if (max_t <= min_t) max_t = min_t + 1;
    return min_t + (int)(lcg_float(&ch->rng) * (float)(max_t - min_t));
}

static void chew_set_params(chew_t *ch, float chew_param) {
    ch->depth = chew_param;
    ch->freq = chew_param;
    ch->variance = chew_param * 0.5f;
}

static void chew_process_block(chew_t *ch, float *left, float *right, int frames) {
    float highFreq = fminf(22000.0f, 0.49f * ch->fs);
    float freqChange = highFreq - 5000.0f;

    /* Process in 64-sample chunks for parameter updates */
    for (int start = 0; start < frames; start += 64) {
        int chunk = frames - start;
        if (chunk > 64) chunk = 64;

        /* Update chew state */
        if (ch->freq == 0.0f) {
            smooth_set(&ch->dropout_mix, 0.0f);
            onepole_set_freq(&ch->filt[0], highFreq, ch->fs);
            onepole_set_freq(&ch->filt[1], highFreq, ch->fs);
        } else if (ch->sample_counter >= ch->samples_until_change) {
            ch->sample_counter = 0;
            ch->is_crinkled = !ch->is_crinkled;

            if (ch->is_crinkled) {
                smooth_set(&ch->dropout_mix, 1.0f);
                float power = (1.0f + 2.0f * lcg_float(&ch->rng)) * ch->depth;
                smooth_set(&ch->dropout_power, 1.0f + power);
                float filterFreq = highFreq - freqChange * ch->depth;
                onepole_set_freq(&ch->filt[0], filterFreq, ch->fs);
                onepole_set_freq(&ch->filt[1], filterFreq, ch->fs);
                ch->samples_until_change = chew_get_wet_time(ch);
            } else {
                smooth_set(&ch->dropout_mix, 0.0f);
                onepole_set_freq(&ch->filt[0], highFreq, ch->fs);
                onepole_set_freq(&ch->filt[1], highFreq, ch->fs);
                ch->samples_until_change = chew_get_dry_time(ch);
            }
        } else if (ch->is_crinkled) {
            float power = (1.0f + 2.0f * lcg_float(&ch->rng)) * ch->depth;
            smooth_set(&ch->dropout_power, 1.0f + power);
            float filterFreq = highFreq - freqChange * ch->depth;
            onepole_set_freq(&ch->filt[0], filterFreq, ch->fs);
            onepole_set_freq(&ch->filt[1], filterFreq, ch->fs);
        }

        /* Process dropout + filter */
        for (int n = 0; n < chunk; n++) {
            int i = start + n;
            float mix = smooth_next(&ch->dropout_mix);
            float pwr = smooth_next(&ch->dropout_power);

            if (mix > 0.001f) {
                /* Dropout: power distortion */
                float sign_l = (left[i] >= 0.0f) ? 1.0f : -1.0f;
                float sign_r = (right[i] >= 0.0f) ? 1.0f : -1.0f;
                float drop_l = powf(fabsf(left[i]), pwr) * sign_l;
                float drop_r = powf(fabsf(right[i]), pwr) * sign_r;
                left[i] = mix * drop_l + (1.0f - mix) * left[i];
                right[i] = mix * drop_r + (1.0f - mix) * right[i];
            }

            left[i] = onepole_process(&ch->filt[0], left[i]);
            right[i] = onepole_process(&ch->filt[1], right[i]);
        }

        ch->sample_counter += chunk;
    }
}

/* ==================================================================
 * Degrade - Noise + lowpass
 * ================================================================== */

typedef struct {
    onepole_t filt[2];
    uint32_t rng;
    float noise_gain;
    float filter_freq;
    float gain_linear;
    float fs;
} degrade_t;

static void degrade_init(degrade_t *d, float fs) {
    onepole_reset(&d->filt[0]);
    onepole_reset(&d->filt[1]);
    onepole_set_freq(&d->filt[0], 20000.0f, fs);
    onepole_set_freq(&d->filt[1], 20000.0f, fs);
    d->rng = 54321;
    d->noise_gain = 0.0f;
    d->filter_freq = 20000.0f;
    d->gain_linear = 1.0f;
    d->fs = fs;
}

static void degrade_set_params(degrade_t *d, float degrade_param) {
    if (degrade_param < 0.01f) {
        d->noise_gain = 0.0f;
        d->gain_linear = 1.0f;
        onepole_set_freq(&d->filt[0], 20000.0f, d->fs);
        onepole_set_freq(&d->filt[1], 20000.0f, d->fs);
        return;
    }

    /* Noise increases with degradation */
    d->noise_gain = 0.5f * degrade_param * degrade_param;

    /* Filter frequency decreases with degradation */
    float freq_hz = 200.0f * powf(20000.0f / 200.0f, 1.0f - degrade_param);
    d->filter_freq = fminf(freq_hz, 0.49f * d->fs);
    onepole_set_freq(&d->filt[0], d->filter_freq, d->fs);
    onepole_set_freq(&d->filt[1], d->filter_freq, d->fs);

    /* Gain loss */
    float gain_db = -24.0f * degrade_param;
    d->gain_linear = powf(10.0f, gain_db / 20.0f);
}

static void degrade_process(degrade_t *d, float *left, float *right, int frames) {
    if (d->noise_gain < 0.001f && d->gain_linear > 0.999f) return;

    for (int n = 0; n < frames; n++) {
        /* Add noise */
        float noise_l = (lcg_float(&d->rng) - 0.5f) * d->noise_gain;
        float noise_r = (lcg_float(&d->rng) - 0.5f) * d->noise_gain;
        left[n] += noise_l;
        right[n] += noise_r;

        /* Lowpass filter */
        left[n] = onepole_process(&d->filt[0], left[n]);
        right[n] = onepole_process(&d->filt[1], right[n]);

        /* Gain reduction */
        left[n] *= d->gain_linear;
        right[n] *= d->gain_linear;
    }
}

/* ==================================================================
 * Loss Filter - Head loss modeled as lowpass
 * Speed -> cutoff frequency mapping based on tape physics
 * ================================================================== */

typedef struct {
    onepole_t filt[2];
    /* Head bump resonance */
    biquad_t bump[2];
    float fs;
} lossfilter_t;

static void loss_init(lossfilter_t *lf, float fs) {
    onepole_reset(&lf->filt[0]);
    onepole_reset(&lf->filt[1]);
    biquad_reset(&lf->bump[0]);
    biquad_reset(&lf->bump[1]);
    lf->fs = fs;
}

static void loss_set_speed(lossfilter_t *lf, float speed_param) {
    /* Speed param 0-1 maps to tape speed 1-50 IPS (inches per second) */
    float speed_ips = 1.0f + 49.0f * speed_param;

    /* Loss filter cutoff: higher speed = higher cutoff */
    /* Based on spacing loss: H = exp(-k*d), cutoff roughly proportional to speed */
    float cutoff = 2000.0f + 18000.0f * powf(speed_param, 0.7f);
    cutoff = fminf(cutoff, 0.49f * lf->fs);

    onepole_set_freq(&lf->filt[0], cutoff, lf->fs);
    onepole_set_freq(&lf->filt[1], cutoff, lf->fs);

    /* Head bump filter: resonance at low frequency proportional to speed */
    float bump_freq = speed_ips * 0.0254f / (1.0e-6f * 500.0f);
    bump_freq = fmaxf(20.0f, fminf(bump_freq, 1000.0f));
    float bump_gain = fmaxf(1.5f * (1000.0f - fabsf(bump_freq - 100.0f)) / 1000.0f, 1.0f);
    float bump_db = 20.0f * log10f(bump_gain);

    /* Peak filter for head bump */
    float w0 = (float)kTwoPi * bump_freq / lf->fs;
    float sinw0 = sinf(w0);
    float cosw0 = cosf(w0);
    float A = powf(10.0f, bump_db / 40.0f);
    float alpha = sinw0 / (2.0f * 2.0f); /* Q = 2 */

    float a0 = 1.0f + alpha / A;
    for (int ch = 0; ch < 2; ch++) {
        lf->bump[ch].b0 = (1.0f + alpha * A) / a0;
        lf->bump[ch].b1 = (-2.0f * cosw0) / a0;
        lf->bump[ch].b2 = (1.0f - alpha * A) / a0;
        lf->bump[ch].a1 = (-2.0f * cosw0) / a0;
        lf->bump[ch].a2 = (1.0f - alpha / A) / a0;
    }
}

static inline void loss_process(lossfilter_t *lf, float *l, float *r) {
    *l = onepole_process(&lf->filt[0], *l);
    *r = onepole_process(&lf->filt[1], *r);
    *l = biquad_process(&lf->bump[0], *l);
    *r = biquad_process(&lf->bump[1], *r);
}

/* ==================================================================
 * Instance - Full CHOWTape instance
 * ================================================================== */

typedef struct {
    char module_dir[256];

    /* User parameters (0-1 normalized) */
    float param_drive;
    float param_saturation;
    float param_bias;
    float param_tone;
    float param_speed;
    float param_flutter;
    float param_chew;
    float param_degrade;
    float param_mix;
    float param_output;

    /* Hysteresis (one per channel) */
    hysteresis_t hyst_l;
    hysteresis_t hyst_r;

    /* Tone */
    biquad_t tone_low[2];
    biquad_t tone_high[2];

    /* Loss filter */
    lossfilter_t loss;

    /* Wow & Flutter */
    wowflutter_t wf;

    /* Chew */
    chew_t chew;

    /* Degrade */
    degrade_t degrade;

    /* DC Blocker (post-hysteresis) */
    dc_blocker_t dc_post[2];

    float fs;
} chowtape_instance_t;

static void update_tone(chowtape_instance_t *inst) {
    /* tone 0-1: 0=dark (-6dB treble, +3dB bass), 0.5=flat, 1=bright (+6dB treble, -3dB bass) */
    float tone = inst->param_tone;
    float treble_db = (tone - 0.5f) * 12.0f;
    float bass_db = -(tone - 0.5f) * 6.0f;
    float trans_freq = 500.0f;

    for (int ch = 0; ch < 2; ch++) {
        biquad_lowshelf(&inst->tone_low[ch], inst->fs, trans_freq, bass_db);
        biquad_highshelf(&inst->tone_high[ch], inst->fs, trans_freq, treble_db);
    }
}

static void update_hysteresis(chowtape_instance_t *inst) {
    hyst_cook(&inst->hyst_l, (double)inst->param_drive, (double)inst->param_bias, (double)inst->param_saturation);
    hyst_cook(&inst->hyst_r, (double)inst->param_drive, (double)inst->param_bias, (double)inst->param_saturation);
}

/* ---- Audio FX API v2 ---- */

typedef struct audio_fx_api_v2 {
    uint32_t api_version;
    void* (*create_instance)(const char *module_dir, const char *config_json);
    void (*destroy_instance)(void *instance);
    void (*process_block)(void *instance, int16_t *audio_inout, int frames);
    void (*set_param)(void *instance, const char *key, const char *val);
    int (*get_param)(void *instance, const char *key, char *buf, int buf_len);
} audio_fx_api_v2_t;

static void* v2_create_instance(const char *module_dir, const char *config_json) {
    chowtape_instance_t *inst = (chowtape_instance_t*)calloc(1, sizeof(chowtape_instance_t));
    if (!inst) return NULL;

    if (module_dir)
        strncpy(inst->module_dir, module_dir, sizeof(inst->module_dir) - 1);

    inst->fs = (float)SAMPLE_RATE;

    /* Default parameters */
    inst->param_drive = 0.5f;
    inst->param_saturation = 0.5f;
    inst->param_bias = 0.5f;
    inst->param_tone = 0.5f;
    inst->param_speed = 0.5f;
    inst->param_flutter = 0.0f;
    inst->param_chew = 0.0f;
    inst->param_degrade = 0.0f;
    inst->param_mix = 1.0f;
    inst->param_output = 0.5f;

    /* Init DSP modules */
    hyst_init(&inst->hyst_l, SAMPLE_RATE);
    hyst_init(&inst->hyst_r, SAMPLE_RATE);
    update_hysteresis(inst);

    for (int ch = 0; ch < 2; ch++) {
        biquad_reset(&inst->tone_low[ch]);
        biquad_reset(&inst->tone_high[ch]);
    }
    update_tone(inst);

    loss_init(&inst->loss, inst->fs);
    loss_set_speed(&inst->loss, inst->param_speed);

    wf_init(&inst->wf, inst->fs);
    wf_set_params(&inst->wf, inst->param_flutter);

    chew_init(&inst->chew, inst->fs);
    chew_set_params(&inst->chew, inst->param_chew);

    degrade_init(&inst->degrade, inst->fs);
    degrade_set_params(&inst->degrade, inst->param_degrade);

    dcblock_init(&inst->dc_post[0], 15.0f, inst->fs);
    dcblock_init(&inst->dc_post[1], 15.0f, inst->fs);

    host_log("[chowtape] instance created");
    return inst;
}

static void v2_destroy_instance(void *instance) {
    if (instance) {
        host_log("[chowtape] instance destroyed");
        free(instance);
    }
}

static inline float soft_clip(float x) {
    if (x > 1.0f) return 1.0f;
    if (x < -1.0f) return -1.0f;
    return 1.5f * x - 0.5f * x * x * x;
}

static void v2_process_block(void *instance, int16_t *audio_inout, int frames) {
    chowtape_instance_t *inst = (chowtape_instance_t*)instance;
    if (!inst) return;

    float left_buf[512], right_buf[512];
    float dry_l[512], dry_r[512];

    float output_gain = inst->param_output * 2.0f;  /* 0-1 -> 0-2x */
    float mix = inst->param_mix;

    for (int chunk = 0; chunk < frames; chunk += 512) {
        int n_frames = frames - chunk;
        if (n_frames > 512) n_frames = 512;

        int16_t *buf = audio_inout + chunk * 2;

        /* Convert int16 to float and save dry signal */
        for (int n = 0; n < n_frames; n++) {
            float l = (float)buf[n * 2] / 32768.0f;
            float r = (float)buf[n * 2 + 1] / 32768.0f;
            dry_l[n] = l;
            dry_r[n] = r;
            left_buf[n] = l;
            right_buf[n] = r;
        }

        /* === Processing chain === */

        /* 1. Tone (pre-EQ, like the original's input tone stage) */
        for (int n = 0; n < n_frames; n++) {
            left_buf[n] = biquad_process(&inst->tone_low[0], left_buf[n]);
            left_buf[n] = biquad_process(&inst->tone_high[0], left_buf[n]);
            right_buf[n] = biquad_process(&inst->tone_low[1], right_buf[n]);
            right_buf[n] = biquad_process(&inst->tone_high[1], right_buf[n]);
        }

        /* 2. Hysteresis (the core tape saturation) */
        for (int n = 0; n < n_frames; n++) {
            /* Scale input for hysteresis (the model expects small values) */
            double in_l = (double)left_buf[n];
            double in_r = (double)right_buf[n];

            double out_l = hyst_process(&inst->hyst_l, in_l);
            double out_r = hyst_process(&inst->hyst_r, in_r);

            left_buf[n] = (float)out_l;
            right_buf[n] = (float)out_r;

            /* DC blocking after hysteresis */
            left_buf[n] = dcblock_process(&inst->dc_post[0], left_buf[n]);
            right_buf[n] = dcblock_process(&inst->dc_post[1], right_buf[n]);
        }

        /* 3. Loss filter (head loss / HF rolloff) */
        for (int n = 0; n < n_frames; n++) {
            loss_process(&inst->loss, &left_buf[n], &right_buf[n]);
        }

        /* 4. Chew (tape damage) */
        if (inst->param_chew > 0.01f) {
            chew_process_block(&inst->chew, left_buf, right_buf, n_frames);
        }

        /* 5. Degrade (noise + filter) */
        if (inst->param_degrade > 0.01f) {
            degrade_process(&inst->degrade, left_buf, right_buf, n_frames);
        }

        /* 6. Wow & Flutter */
        if (inst->param_flutter > 0.01f) {
            for (int n = 0; n < n_frames; n++) {
                wf_process_sample(&inst->wf, &left_buf[n], &right_buf[n]);
            }
        }

        /* 7. Mix dry/wet + output gain + convert back to int16 */
        for (int n = 0; n < n_frames; n++) {
            float l = mix * left_buf[n] + (1.0f - mix) * dry_l[n];
            float r = mix * right_buf[n] + (1.0f - mix) * dry_r[n];

            l *= output_gain;
            r *= output_gain;

            /* Soft clip */
            l = soft_clip(l);
            r = soft_clip(r);

            /* Convert to int16 */
            int32_t il = (int32_t)(l * 32767.0f);
            int32_t ir = (int32_t)(r * 32767.0f);
            if (il > 32767) il = 32767;
            if (il < -32767) il = -32767;
            if (ir > 32767) ir = 32767;
            if (ir < -32767) ir = -32767;

            buf[n * 2] = (int16_t)il;
            buf[n * 2 + 1] = (int16_t)ir;
        }
    }
}

static void v2_set_param(void *instance, const char *key, const char *val) {
    chowtape_instance_t *inst = (chowtape_instance_t*)instance;
    if (!inst || !key || !val) return;

    float fval = (float)atof(val);
    if (fval < 0.0f) fval = 0.0f;
    if (fval > 1.0f) fval = 1.0f;

    if (strcmp(key, "drive") == 0) {
        inst->param_drive = fval;
        update_hysteresis(inst);
    } else if (strcmp(key, "saturation") == 0) {
        inst->param_saturation = fval;
        update_hysteresis(inst);
    } else if (strcmp(key, "bias") == 0) {
        inst->param_bias = fval;
        update_hysteresis(inst);
    } else if (strcmp(key, "tone") == 0) {
        inst->param_tone = fval;
        update_tone(inst);
    } else if (strcmp(key, "speed") == 0) {
        inst->param_speed = fval;
        loss_set_speed(&inst->loss, fval);
    } else if (strcmp(key, "flutter") == 0) {
        inst->param_flutter = fval;
        wf_set_params(&inst->wf, fval);
    } else if (strcmp(key, "chew") == 0) {
        inst->param_chew = fval;
        chew_set_params(&inst->chew, fval);
    } else if (strcmp(key, "degrade") == 0) {
        inst->param_degrade = fval;
        degrade_set_params(&inst->degrade, fval);
    } else if (strcmp(key, "mix") == 0) {
        inst->param_mix = fval;
    } else if (strcmp(key, "output") == 0) {
        inst->param_output = fval;
    } else if (strcmp(key, "state") == 0) {
        /* Restore state from JSON - parse key:value pairs */
        const char *p = val;
        while (*p) {
            /* Find key */
            const char *kstart = strchr(p, '"');
            if (!kstart) break;
            kstart++;
            const char *kend = strchr(kstart, '"');
            if (!kend) break;

            /* Find value */
            const char *vstart = strchr(kend + 1, ':');
            if (!vstart) break;
            vstart++;
            while (*vstart == ' ' || *vstart == '"') vstart++;

            char kbuf[64], vbuf[64];
            int klen = (int)(kend - kstart);
            if (klen >= 64) klen = 63;
            memcpy(kbuf, kstart, klen);
            kbuf[klen] = '\0';

            /* Extract value until comma, }, or quote */
            int vi = 0;
            while (vstart[vi] && vstart[vi] != ',' && vstart[vi] != '}' && vstart[vi] != '"' && vi < 63) vi++;
            memcpy(vbuf, vstart, vi);
            vbuf[vi] = '\0';

            if (strcmp(kbuf, "state") != 0) {
                v2_set_param(instance, kbuf, vbuf);
            }

            p = vstart + vi;
        }
    }
}

static int v2_get_param(void *instance, const char *key, char *buf, int buf_len) {
    chowtape_instance_t *inst = (chowtape_instance_t*)instance;
    if (!inst || !key || !buf || buf_len < 1) return -1;

    if (strcmp(key, "drive") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_drive);
    if (strcmp(key, "saturation") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_saturation);
    if (strcmp(key, "bias") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_bias);
    if (strcmp(key, "tone") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_tone);
    if (strcmp(key, "speed") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_speed);
    if (strcmp(key, "flutter") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_flutter);
    if (strcmp(key, "chew") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_chew);
    if (strcmp(key, "degrade") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_degrade);
    if (strcmp(key, "mix") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_mix);
    if (strcmp(key, "output") == 0)
        return snprintf(buf, buf_len, "%.2f", inst->param_output);

    if (strcmp(key, "name") == 0)
        return snprintf(buf, buf_len, "CHOWTape");

    if (strcmp(key, "state") == 0) {
        return snprintf(buf, buf_len,
            "{\"drive\":%.2f,\"saturation\":%.2f,\"bias\":%.2f,\"tone\":%.2f,"
            "\"speed\":%.2f,\"flutter\":%.2f,\"chew\":%.2f,\"degrade\":%.2f,"
            "\"mix\":%.2f,\"output\":%.2f}",
            inst->param_drive, inst->param_saturation, inst->param_bias,
            inst->param_tone, inst->param_speed, inst->param_flutter,
            inst->param_chew, inst->param_degrade, inst->param_mix,
            inst->param_output);
    }

    if (strcmp(key, "chain_params") == 0) {
        return snprintf(buf, buf_len,
            "[{\"key\":\"drive\",\"label\":\"Drive\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0.5,\"step\":0.01},"
            "{\"key\":\"saturation\",\"label\":\"Saturation\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0.5,\"step\":0.01},"
            "{\"key\":\"bias\",\"label\":\"Bias\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0.5,\"step\":0.01},"
            "{\"key\":\"tone\",\"label\":\"Tone\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0.5,\"step\":0.01},"
            "{\"key\":\"speed\",\"label\":\"Speed\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0.5,\"step\":0.01},"
            "{\"key\":\"flutter\",\"label\":\"Flutter\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0,\"step\":0.01},"
            "{\"key\":\"chew\",\"label\":\"Chew\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0,\"step\":0.01},"
            "{\"key\":\"degrade\",\"label\":\"Degrade\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0,\"step\":0.01},"
            "{\"key\":\"mix\",\"label\":\"Mix\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":1,\"step\":0.01},"
            "{\"key\":\"output\",\"label\":\"Output\",\"type\":\"float\",\"min\":0,\"max\":1,\"default\":0.5,\"step\":0.01}]");
    }

    if (strcmp(key, "ui_hierarchy") == 0) {
        return snprintf(buf, buf_len,
            "{\"levels\":{\"root\":{\"name\":\"CHOWTape\","
            "\"knobs\":[\"drive\",\"saturation\",\"bias\",\"tone\",\"speed\",\"flutter\",\"chew\",\"mix\"]}}}");
    }

    return -1;
}

/* ---- Module entry point ---- */

static audio_fx_api_v2_t g_api;

__attribute__((visibility("default")))
audio_fx_api_v2_t* move_audio_fx_init_v2(const host_api_v1_t *host) {
    g_host = host;

    g_api.api_version = AUDIO_FX_API_VERSION_2;
    g_api.create_instance = v2_create_instance;
    g_api.destroy_instance = v2_destroy_instance;
    g_api.process_block = v2_process_block;
    g_api.set_param = v2_set_param;
    g_api.get_param = v2_get_param;

    if (host && host->log)
        host->log("[chowtape] CHOWTape v0.1.0 loaded (GPL-3.0, based on AnalogTapeModel by Jatin Chowdhury)");

    return &g_api;
}

/*
    Splat - _splat.h

    Copyright (C) 2015, 2016, 2017
    Guillaume Tucker <guillaume@mangoz.org>

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU Lesser General Public License as published by
    the Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
    License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _SPLAT_H
#define _SPLAT_H 1

#include <Python.h>
#include <math.h>

/* Enable for speed using SIMD and 32-bit samples instead of 64-bit */
#ifdef SPLAT_FAST
#if defined(__ARM_NEON__) || defined(__ARM_NEON)
#include <arm_neon.h>
#define SPLAT_NEON
typedef float32x4_t sf_float_t;
typedef uint32x4_t sf_uint_t;
typedef uint32x4_t sf_mask_t;
#define sf_set(x) vdupq_n_f32((x))
#define sf_zero() vdupq_n_f32(0.0)
#define sf_add(x, y) vaddq_f32((x), (y))
#define sf_mul(x, y) vmulq_f32((x), (y))
#define sf_and(x, y) (sf_float_t)vandq_u32((sf_uint_t)(x), (sf_uint_t)(y))
#define sf_lt(x, y) vcltq_f32((x), (y))
#elif defined(__SSE__)
#include <x86intrin.h>
#define SPLAT_SSE
typedef __m128 sf_float_t;
typedef __m128i sf_uint_t;
typedef __m128 sf_mask_t;
#define sf_set(x) _mm_set1_ps((x))
#define sf_zero() _mm_setzero_ps()
#define sf_add(x, y) _mm_add_ps((x), (y))
#define sf_mul(x, y) _mm_mul_ps((x), (y))
#define sf_and(x, y) _mm_and_ps((x), (y))
#define sf_lt(x, y) _mm_cmplt_ps((x), (y))
#endif
#endif

/* Sample type */
#ifdef SPLAT_FAST
typedef float sample_t;
#define SPLAT_NATIVE_SAMPLE_TYPE SPLAT_FLOAT_32
#else
typedef double sample_t;
#define SPLAT_NATIVE_SAMPLE_TYPE SPLAT_FLOAT_64
#endif

/* Maximum number of channels */
#define SPLAT_MAX_CHANNELS 16

/* Linear and dB conversions */
#define lin2dB(level) (20 * log10(level))
#define dB2lin(dB) (pow10((dB) / 20))

/* Round to the next multiple of 4, used in fast mode for 4x float vectors */
#ifdef SPLAT_FAST
#define splat_round4(_len) ((_len) % 4 ? ((((_len) / 4) + 1) * 4) : (_len))
#define splat_mask4(_x) (_x & 0xFFFFFFFFFFFFFFFC)
#endif

#ifndef min
#define min(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif

#ifndef max
#define max(_a, _b) (((_a) > (_b)) ? (_a) : (_b))
#endif

#ifndef minmax
#define minmax(_val, _min, _max) \
	((_val) < (_min) ? (_min) : ((_val) > (_max) ? (_max) : (_val)))
#endif

#ifndef ARRAY_SIZE
#define ARRAY_SIZE(_array) (sizeof(_array) / sizeof(_array[0]))
#endif

/* Page utilities */
extern size_t splat_page_size;
extern const void *splat_zero_page;

/* Convert any number type to a double or return -1 */
extern int splat_obj2double(PyObject *obj, double *out);

/* Levels (or gains) */
struct splat_levels {
	unsigned n;
	PyObject *obj[SPLAT_MAX_CHANNELS];
	double fl[SPLAT_MAX_CHANNELS]; /* levels converted to linear scale */
	int all_floats;
#ifdef SPLAT_FAST
	sf_float_t flq[SPLAT_MAX_CHANNELS];
#endif
};

/* Fast sine function interpolation table */
struct splat_sine_poly {
	float coef[4];
};

extern const struct splat_sine_poly *splat_sine_table;
extern const size_t splat_sine_table_len;
extern const size_t splat_sine_table_mask;

#ifdef SPLAT_FAST
#define SPLAT_QUAD(_x) { (_x), (_x), (_x), (_x) }
extern sf_float_t splat_fast_sine_step;
extern sf_mask_t splat_fast_sine_mask;
extern const sf_float_t splat_fast_inc;
#endif

#if defined(SPLAT_SSE)
extern const __m128 splat_sse_zero;
extern const __m128 splat_sse_one;
extern const __m128 splat_sse_two;
extern const __m128 splat_sse_pi;
extern const __m128 splat_sse_pi2;
extern const __m128 splat_sse_inc;
#endif

/* ----------------------------------------------------------------------------
 * mmap
 */

struct splat_mmap {
	char *path;
	int fd;
	void *ptr;
	size_t size;
	int persist;
};

extern int splat_mmap_init(struct splat_mmap *m, const char *path);
extern void splat_mmap_free(struct splat_mmap *m);
extern int splat_mmap_remap(struct splat_mmap *m, size_t length);
extern int splat_mmap_set_temp_px(const char *px);
extern const char *splat_mmap_get_temp_px(void);

/* ----------------------------------------------------------------------------
 * Fragment
 */

struct splat_channel {
	sample_t *data;
	size_t length;
	struct splat_mmap mmap;
	void (*free)(struct splat_channel *chan);
	int (*resize)(struct splat_channel *chan, size_t length);
};

struct splat_fragment {
	size_t length;
	unsigned rate;
	struct splat_channel channels[SPLAT_MAX_CHANNELS];
	unsigned n_channels;
	char *name;
	int uses_mmap;
	int temp_mmap;
};

struct splat_peak {
	double avg;
	double max;
	double min;
	double peak;
};

struct splat_levels;

extern struct splat_fragment *splat_frag_from_obj(PyObject *obj);
extern int splat_frag_init(struct splat_fragment *frag, unsigned n_channels,
			   unsigned rate, size_t length, const char *name);
extern int splat_frag_init_mmap(struct splat_fragment *frag,
				unsigned n_channels, unsigned rate,
				size_t length, const char *name,
				const char *new_path, const char **open_paths);
extern void splat_frag_free(struct splat_fragment *frag);
extern int splat_frag_set_name(struct splat_fragment *frag, const char *name);
extern int splat_frag_resize(struct splat_fragment *frag, size_t length);
#define splat_frag_grow(_frag, _length)		\
	(((_length) <= (_frag)->length) ? 0 :	\
	 splat_frag_resize((_frag), (_length)))
extern int splat_frag_mix(struct splat_fragment *frag,
			  const struct splat_fragment *incoming,
			  const struct splat_levels *levels, size_t length,
			  double offset, double skip, int zero_dB);
extern int splat_frag_sample_number(size_t *val, long min_val,
				    long max_val, PyObject *obj);
extern void splat_frag_get_peak(const struct splat_fragment *frag,
				struct splat_peak *chan_peak,
				struct splat_peak *frag_peak, int do_avg);
extern void splat_frag_normalize(struct splat_fragment *frag, double level_dB,
				 int do_zero);
extern int splat_frag_amp(struct splat_fragment *frag,
			  struct splat_levels *gains);
extern void splat_frag_lin2dB(struct splat_fragment *frag);
extern void splat_frag_dB2lin(struct splat_fragment *frag);
extern int splat_frag_offset(struct splat_fragment *frag, PyObject *offset_obj,
			     double start);
extern int splat_frag_resample(struct splat_fragment *frag, unsigned rate,
			       PyObject *ratio);

/* ----------------------------------------------------------------------------
 * Signal & vector
 */

#define SPLAT_VECTOR_BITS 8
#define SPLAT_VECTOR_LEN (1 << SPLAT_VECTOR_BITS)

struct splat_signal;

struct splat_vector {
	sample_t data[SPLAT_VECTOR_LEN];
	PyObject *obj;
	int (*signal)(struct splat_signal *s, struct splat_vector *v);
};

enum signal_ret {
	SPLAT_SIGNAL_CONTINUE = 0,
	SPLAT_SIGNAL_ERROR,
	SPLAT_SIGNAL_STOP,
};

struct splat_signal {
	enum signal_ret stat;
	size_t origin;
	size_t length;
	size_t n_vectors;
	struct splat_vector *vectors;
	unsigned rate;
	PyObject *py_float;
	PyObject *py_args;
	size_t cur;
	size_t end;
	size_t len;
};

extern int splat_signal_init(struct splat_signal *s, size_t length,
			     size_t origin, PyObject **signals,
			     size_t n_signals, unsigned rate);
extern void splat_signal_free(struct splat_signal *s);
extern int splat_signal_next(struct splat_signal *s);
extern ssize_t splat_signal_get(struct splat_signal *s, size_t n);
extern PyObject *splat_signal_tuple(struct splat_signal *s, size_t offset);

/* ----------------------------------------------------------------------------
 * Spline
 */

struct splat_spline {
	PyObject *pols; /* List of tuples with (x0, x1, coeffs) */
	double k0;
	double start;
	double end;
	int db;
};

extern struct splat_spline *splat_spline_from_obj(PyObject *obj);
extern double splat_spline_tuple_value(PyObject *poly, double x, int db);
extern PyObject *splat_spline_find_poly(PyObject *spline, double x,
					double *end);

/* ----------------------------------------------------------------------------
 * Sources
 */

struct splat_overtone {
	PyObject *ratio;
	double fl_ratio;
	PyObject *phase;
	double fl_phase;
	struct splat_levels levels;
#ifdef SPLAT_FAST
	sf_float_t fl_ratioq;
	sf_float_t fl_phaseq;
#endif
};

extern void splat_sine_floats(struct splat_fragment *frag,
			      const double *levels, double freq, double phase);
extern int splat_sine_signals(struct splat_fragment *frag, PyObject **levels,
			      PyObject *freq, PyObject *phase, double origin);
extern void splat_square_floats(struct splat_fragment *frag,
				const double *fl_pos, double freq,
				double phase, double ratio);
extern int splat_square_signals(struct splat_fragment *frag, PyObject **levels,
				PyObject *freq, PyObject *phase,
				PyObject *ratio, double origin);
extern void splat_triangle_floats(struct splat_fragment *frag,
				  const double *lvls, double freq,
				  double phase, double ratio);
extern int splat_triangle_signals(struct splat_fragment *frag,
				  PyObject **levels, PyObject *freq,
				  PyObject *phase, PyObject *ratio,
				  double origin);
extern void splat_overtones_float(struct splat_fragment *frag,
				  const double *levels, double freq,
				  double phase,
				  struct splat_overtone *overtones,
				  Py_ssize_t n);
extern int splat_overtones_mixed(struct splat_fragment *frag, PyObject **levels,
				 PyObject *freq, PyObject *phase,
				 struct splat_overtone *overtones,
				 Py_ssize_t n, double origin);
extern int splat_overtones_signal(struct splat_fragment *frag,
				  PyObject **levels, PyObject *freq,
				  PyObject *phase,
				  struct splat_overtone *overtones,
				  Py_ssize_t n, double origin);

/* ----------------------------------------------------------------------------
 * Filters
 */

struct splat_delay {
	size_t time;
	double gain;
};

extern void splat_filter_dec_envelope(struct splat_fragment *frag,
			       double k, double p);
extern void splat_filter_reverse(struct splat_fragment *frag);
extern void splat_filter_reverb(struct splat_fragment *frag,
				struct splat_delay **delays,
				size_t n_delays, size_t max_index);

#endif /* _SPLAT_H */

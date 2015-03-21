/*
    Splat - source.c

    Copyright (C) 2015
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

#include "_splat.h"

#if defined(SPLAT_NEON)
static float32x4_t splat_sine_neon(float32x4_t x);
#elif defined(SPLAT_SSE)
static __m128 splat_sine_sse(__m128 x);
#endif

/* -- sine source -- */

enum {
	SIG_SINE_FREQ = 0,
	SIG_SINE_PHASE,
	SIG_SINE_AMP,
} splat_sine_signal;

#if defined(SPLAT_NEON)
void splat_sine_floats(struct splat_fragment *frag, const double *levels,
		       double freq, double phase)
{
	const float32x4_t phrate = vdupq_n_f32(phase * frag->rate);
	const float32x4_t k = vdupq_n_f32(2.0 * M_PI * freq / frag->rate);
	float32x4_t levelsq[SPLAT_MAX_CHANNELS];
	float32x4_t *out[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	for (c = 0; c < frag->n_channels; ++c) {
		levelsq[c] = vdupq_n_f32(levels[c]);
		out[c] = (float32x4_t *)frag->data[c];
	}

	for (i = 0; i < frag->length; i += 4) {
		uint32x4_t iq;
		float32x4_t x;
		float32x4_t y;

		/* x = 2 * M_PI * freq * (phrate + i) / frag->rate) */
		iq = vdupq_n_u32(i);
		iq = vaddq_u32(iq, splat_neon_inc);
		x = vcvtq_f32_u32(iq);
		x = vaddq_f32(x, phrate);
		x = vmulq_f32(x, k);
		y = splat_sine_neon(x);

		for (c = 0; c < frag->n_channels; ++c)
			*out[c]++ = vmulq_f32(y, levelsq[c]);
	}
}
#elif defined(SPLAT_SSE)
void splat_sine_floats(struct splat_fragment *frag, const double *levels,
		       double freq, double phase)
{
	const __m128 k = _mm_set1_ps(2.0 * M_PI * freq / frag->rate);
	const __m128 phrate = _mm_set1_ps(phase * frag->rate);
	__m128 levelsq[SPLAT_MAX_CHANNELS];
	__m128 *out[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	for (c = 0; c < frag->n_channels; ++c) {
		levelsq[c] = _mm_set1_ps(levels[c]);
		out[c] = (__m128 *)frag->data[c];
	}

	for (i = 0; i < frag->length; i += 4) {
		__m128 x;
		__m128 y;

		x = _mm_set1_ps((float)i);
		x = _mm_add_ps(x, splat_sse_inc);
		x = _mm_add_ps(x, phrate);
		x = _mm_mul_ps(x, k);
		y = splat_sine_sse(x);

		for (c = 0; c < frag->n_channels; ++c)
			*out[c]++ = _mm_mul_ps(y, levelsq[c]);
	}
}
#else
void splat_sine_floats(struct splat_fragment *frag, const double *levels,
		       double freq, double phase)
{
	const double k = 2 * M_PI * freq / frag->rate;
	const long ph = phase * frag->rate;
	size_t i;

	for (i = 0; i < frag->length; ++i) {
		const double s = sin(k * (i + ph));
		unsigned c;

		for (c = 0; c < frag->n_channels; ++c)
			frag->data[c][i] = s * levels[c];
	}
}
#endif

#if defined(SPLAT_SSE)
static void _splat_sine_signals(struct splat_fragment *frag,
				struct splat_signal *sig, double origin)
{
	const __m128 k = _mm_set1_ps(2.0 * M_PI);
	const __m128 rateq = _mm_set1_ps(frag->rate);
	const __m128 originq = _mm_set1_ps(origin);
	__m128 *out[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i = 0;

	for (c = 0; c < frag->n_channels; ++c)
		out[c] = (__m128 *)frag->data[c];

	while (splat_signal_next(sig) == SPLAT_SIGNAL_CONTINUE) {
		const __m128 *fq = (__m128 *)sig->vectors[SIG_SINE_FREQ].data;
		const __m128 *phq = (__m128 *)sig->vectors[SIG_SINE_PHASE].data;
		const __m128 *aq[SPLAT_MAX_CHANNELS];
		size_t j;

		for (c = 0; c < frag->n_channels; ++c)
			aq[c] = (__m128 *)sig->vectors[SIG_SINE_AMP + c].data;

		for (j = 0; j < sig->len; i += 4, j += 4) {
			__m128 x;
			__m128 f;
			__m128 y;

			x = _mm_set1_ps((float)i);
			x = _mm_add_ps(x, splat_sse_inc);
			x = _mm_div_ps(x, rateq);
			x = _mm_add_ps(x, *phq++);
			x = _mm_add_ps(x, originq);
			f = _mm_mul_ps(*fq++, k);
			x = _mm_mul_ps(x, f);
			y = splat_sine_sse(x);

			for (c = 0; c < frag->n_channels; ++c)
				*out[c]++ = _mm_mul_ps(y, *aq[c]++);
		}
	}
}
#else
static void _splat_sine_signals(struct splat_fragment *frag,
				struct splat_signal *sig, double origin)
{
	static const double k = 2.0 * M_PI;
	size_t i = 0;

	while (splat_signal_next(sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig->len; ++i, ++j) {
			const double f = sig->vectors[SIG_SINE_FREQ].data[j];
			const double ph = sig->vectors[SIG_SINE_PHASE].data[j];
			const double t = ph + origin + (double)i / frag->rate;
			const double s = sin(k * f * t);
			unsigned c;

			for (c = 0; c < frag->n_channels; ++c) {
				const double a =
					sig->vectors[SIG_SINE_AMP + c].data[j];

				frag->data[c][i] = s * a;
			}
		}
	}
}
#endif

int splat_sine_signals(struct splat_fragment *frag, PyObject **levels,
		       PyObject *freq, PyObject *phase, double origin)
{
	struct splat_signal sig;
	PyObject *signals[SIG_SINE_AMP + SPLAT_MAX_CHANNELS];
	unsigned c;

	signals[SIG_SINE_FREQ] = freq;
	signals[SIG_SINE_PHASE] = phase;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_SINE_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, (SIG_SINE_AMP + frag->n_channels),
			      frag->rate))
		return -1;

	_splat_sine_signals(frag, &sig, origin);

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

/* -- square source -- */

void splat_square_floats(struct splat_fragment *frag, const double *fl_pos,
			 double freq, double phase, double ratio)
{
	const double k = freq / frag->rate;
	double fl_neg[SPLAT_MAX_CHANNELS];
	unsigned c;
	Py_ssize_t i;

	ratio = min(ratio, 1.0);
	ratio = max(ratio, 0.0);

	for (c = 0; c < frag->n_channels; ++c)
		fl_neg[c] = -fl_pos[c];

	for (i = 0; i < frag->length; ++i) {
		double n_periods;
		const double t_rel = modf(((i * k) + phase), &n_periods);
		const double *levels;

		if (t_rel < ratio)
			levels = fl_pos;
		else
			levels = fl_neg;

		for (c = 0; c < frag->n_channels; ++c)
			frag->data[c][i] = levels[c];
	}
}

int splat_square_signals(struct splat_fragment *frag, PyObject **levels,
			 PyObject *freq, PyObject *phase,
			 PyObject *ratio, double origin)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_RATIO,
		SIG_AMP,
	};
	struct splat_signal sig;
	PyObject *signals[SIG_AMP + SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	signals[SIG_FREQ] = freq;
	signals[SIG_PHASE] = phase;
	signals[SIG_RATIO] = ratio;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, (SIG_AMP + frag->n_channels),
			      frag->rate))
		return -1;

	i = 0;

	while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig.len; ++i, ++j) {
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double t = ph + origin + (double)i / frag->rate;
			double n_periods;
			const double t_rel = modf((f * t), &n_periods);
			double ratio = sig.vectors[SIG_RATIO].data[j];
			double s;

			ratio = min(ratio, 1.0);
			ratio = max(ratio, 0.0);
			s = (t_rel < ratio) ? 1.0 : -1.0;

			for (c = 0; c < frag->n_channels; ++c) {
				const double a =
					sig.vectors[SIG_AMP + c].data[j];

				frag->data[c][i] = a * s;
			}
		}
	}

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

/* -- triangle source -- */

void splat_triangle_floats(struct splat_fragment *frag, const double *lvls,
			   double freq, double phase, double ratio)
{
	const double k = freq / frag->rate;
	double a1[SPLAT_MAX_CHANNELS], b1[SPLAT_MAX_CHANNELS];
	double a2[SPLAT_MAX_CHANNELS], b2[SPLAT_MAX_CHANNELS];
	const double *a, *b;
	unsigned c;
	Py_ssize_t i;

	ratio = min(ratio, 1.0);
	ratio = max(ratio, 0.0);

	for (c = 0; c < frag->n_channels; ++c) {
		const double llin = lvls[c];

		a1[c] = 2 * llin / ratio;
		b1[c] = -llin;
		a2[c] = -2 * llin / (1 - ratio);
		b2[c] = llin - (a2[c] * ratio);
	}

	for (i = 0; i < frag->length; ++i) {
		double n_periods;
		const double t_rel = modf(((i * k) + phase), &n_periods);

		if (t_rel < ratio) {
			a = a1;
			b = b1;
		} else {
			a = a2;
			b = b2;
		}

		for (c = 0; c < frag->n_channels; ++c)
			frag->data[c][i] = (a[c] * t_rel) + b[c];
	}
}

int splat_triangle_signals(struct splat_fragment *frag, PyObject **levels,
			   PyObject *freq, PyObject *phase, PyObject *ratio,
			   double origin)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_RATIO,
		SIG_AMP,
	};
	struct splat_signal sig;
	PyObject *signals[SIG_AMP + SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	signals[SIG_FREQ] = freq;
	signals[SIG_PHASE] = phase;
	signals[SIG_RATIO] = ratio;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, (SIG_AMP + frag->n_channels),
			      frag->rate))
		return -1;

	i = 0;

	while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig.len; ++i, ++j) {
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double t = ph + origin + (double)i / frag->rate;
			double n_periods;
			const double t_rel = modf((f * t), &n_periods);
			double ratio = sig.vectors[SIG_RATIO].data[j];

			ratio = min(ratio, 1.0);
			ratio = max(ratio, 0.0);

			for (c = 0; c < frag->n_channels; ++c) {
				const double l =
					sig.vectors[SIG_AMP + c].data[j];
				double a, b;

				if (t_rel < ratio) {
					a = 2 * l / ratio;
					b = -l;
				} else {
					a = -2 * l / (1 - ratio);
					b = l - (a * ratio);
				}

				frag->data[c][i] = (a * t_rel) + b;
			}
		}
	}

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

/* -- overtones source -- */

#if defined(SPLAT_NEON)
static void _splat_overtones_float(struct splat_fragment *frag,
				   double freq, double phase,
				   const struct splat_overtone *overtones,
				   const struct splat_overtone *ot_end)
{
	const float32x4_t k = vdupq_n_f32(2 * M_PI * freq);
	const float32x4_t ph = vdupq_n_f32(phase);
	const float32x4_t rate_inv = vdupq_n_f32(1.0 / frag->rate);
	const float32x4_t max_ratio = vdupq_n_f32((frag->rate / 2.0) / freq);
	float32x4_t *outq[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	for (c = 0; c < frag->n_channels; ++c)
		outq[c] = (float32x4_t *)frag->data[c];

	/* x = 2 * M_PI * freq * ratio * (ph + i/rate) */
	for (i = 0; i < frag->length; i += 4) {
		const struct splat_overtone *ot;
		uint32x4_t iq;
		float32x4_t x0;

		iq = vdupq_n_u32(i);
		iq = vaddq_u32(iq, splat_neon_inc);
		x0 = vcvtq_f32_u32(iq);
		x0 = vmulq_f32(x0, rate_inv);
		x0 = vaddq_f32(x0, ph);

		for (ot = overtones; ot != ot_end; ++ot) {
			const float32x4_t *lvlq = ot->levels.flq;
			float32x4_t xr;
			float32x4_t yr;
			uint32x4_t clip;

			clip = vcltq_f32(ot->fl_ratioq, max_ratio);
			xr = vaddq_f32(x0, ot->fl_phaseq);
			xr = vmulq_f32(xr, k);
			xr = vmulq_f32(xr, ot->fl_ratioq);
			yr = splat_sine_neon(xr);
			yr = (float32x4_t)vandq_u32((uint32x4_t)yr, clip);

			for (c = 0; c < frag->n_channels; ++c) {
				const float32x4_t y = vmulq_f32(yr, lvlq[c]);
				*outq[c] = vaddq_f32(*outq[c], y);
			}
		}

		for (c = 0; c < frag->n_channels; ++c)
			outq[c]++;
	}
}
#elif defined(SPLAT_SSE)
static void _splat_overtones_float(struct splat_fragment *frag,
				   double freq, double phase,
				   const struct splat_overtone *overtones,
				   const struct splat_overtone *ot_end)
{
	const __m128 k = _mm_set1_ps(2.0 * M_PI * freq);
	const __m128 ph = _mm_set1_ps(phase);
	const __m128 rateq = _mm_set1_ps((float)frag->rate);
	const __m128 max_ratio = _mm_set1_ps((frag->rate / 2.0) / freq);
	__m128 *outq[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	for (c = 0; c < frag->n_channels; ++c)
		outq[c] = (__m128 *)frag->data[c];

	/* x = (2 * M_PI * freq * ratio * (ph + i / frag->rate) */
	for (i = 0; i < frag->length; i += 4) {
		const struct splat_overtone *ot;
		__m128 x0;

		x0 = _mm_set1_ps((float)i);
		x0 = _mm_add_ps(x0, splat_sse_inc);
		x0 = _mm_div_ps(x0, rateq);
		x0 = _mm_add_ps(x0, ph);

		for (ot = overtones; ot != ot_end; ++ot) {
			const __m128 *lvlq = ot->levels.flq;
			__m128 clip;
			__m128 x;
			__m128 yr;

			clip = _mm_cmplt_ps(ot->fl_ratioq, max_ratio);
			x = _mm_add_ps(x0, ot->fl_phaseq);
			x = _mm_mul_ps(x, ot->fl_ratioq);
			x = _mm_mul_ps(x, k);
			yr = splat_sine_sse(x);

			for (c = 0; c < frag->n_channels; ++c) {
				const __m128 g = _mm_and_ps(clip, lvlq[c]);
				const __m128 y = _mm_mul_ps(yr, g);

				*outq[c] = _mm_add_ps(*outq[c], y);
			}
		}

		for (c = 0; c < frag->n_channels; ++c)
			outq[c]++;
	}
}
#else
static void _splat_overtones_float(struct splat_fragment *frag,
				   double freq, double phase,
				   const struct splat_overtone *overtones,
				   const struct splat_overtone *ot_end)
{
	const double k = 2 * M_PI * freq;
	size_t i;

	for (i = 0; i < frag->length; ++i) {
		const double t = phase + (double)i / frag->rate;
		const struct splat_overtone *ot;

		for (ot = overtones; ot != ot_end; ++ot) {
			const double s =
				sin(k * ot->fl_ratio * (t + ot->fl_phase));
			unsigned c;

			for (c = 0; c < frag->n_channels; ++c)
				frag->data[c][i] += s * ot->levels.fl[c];
		}
	}
}
#endif

void splat_overtones_float(struct splat_fragment *frag, const double *levels,
			   double freq, double phase,
			   struct splat_overtone *overtones, Py_ssize_t n)
{
	const double max_ratio = (frag->rate / freq) / 2;
	struct splat_overtone *ot;
	const struct splat_overtone * const ot_end = &overtones[n];
	unsigned c;

	/* Silence harmonics above (rate / 2) to avoid spectrum overlap
	   and multiply each overtone levels with global levels. */
	for (ot = overtones; ot != ot_end; ++ot) {
		if (ot->fl_ratio >= max_ratio) {
			for (c = 0; c < frag->n_channels; ++c)
				ot->levels.fl[c] = 0.0f;
		} else {
			for (c = 0; c < frag->n_channels; ++c)
				ot->levels.fl[c] *= levels[c];
		}
	}

	_splat_overtones_float(frag, freq, phase, overtones, ot_end);
}

int splat_overtones_mixed(struct splat_fragment *frag, PyObject **levels,
			  PyObject *freq, PyObject *phase,
			  struct splat_overtone *overtones, Py_ssize_t n,
			  double origin)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_AMP,
	};
	const double k = 2 * M_PI;
	const double half_rate = frag->rate / 2;
	PyObject *signals[SIG_AMP + SPLAT_MAX_CHANNELS];
	struct splat_signal sig;
	struct splat_overtone *ot;
	const struct splat_overtone *ot_end = &overtones[n];
	unsigned c;
	size_t i;

	signals[SIG_FREQ] = freq;
	signals[SIG_PHASE] = phase;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, (SIG_AMP + frag->n_channels),
			      frag->rate))
		return -1;

	i = 0;

	while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig.len; ++i, ++j) {
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double max_ratio = half_rate / f;
			const double m = k * f;
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double t = ph + origin + (double)i / frag->rate;

			for (ot = overtones; ot != ot_end; ++ot) {
				double s;

				if (ot->fl_ratio >= max_ratio)
					continue;

				s = sin(m * ot->fl_ratio * (t + ot->fl_phase));

				for (c = 0; c < frag->n_channels; ++c) {
					double x;

					x = sig.vectors[SIG_AMP + c].data[j];
					x *= ot->levels.fl[c];
					frag->data[c][i] += s * x;
				}
			}
		}
	}

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

int splat_overtones_signal(struct splat_fragment *frag, PyObject **levels,
			   PyObject *freq, PyObject *phase,
			   struct splat_overtone *overtones, Py_ssize_t n,
			   double origin)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE = 1,
		SIG_AMP = 2,
		SIG_OT = 2 + SPLAT_MAX_CHANNELS,
	};
	const double k = 2 * M_PI;
	const double half_rate = frag->rate / 2;
	PyObject **signals;
	struct splat_signal sig;
	struct splat_overtone *ot;
	const struct splat_overtone *ot_end = &overtones[n];
	static const size_t sig_freq = 0;
	static const size_t sig_phase = 1;
	static const size_t sig_amp = 2;
	const size_t sig_ot = sig_amp + frag->n_channels;
	PyObject **sig_ot_it;
	/* for each overtone: ratio, phase and levels */
	const size_t sig_n = sig_ot + (n * (2 + frag->n_channels));
	unsigned c;
	size_t i;

	signals = PyMem_Malloc(sig_n * sizeof(PyObject *));

	if (signals == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	signals[sig_freq] = freq;
	signals[sig_phase] = phase;

	for (c = 0; c < frag->n_channels; ++c)
		signals[sig_amp + c] = levels[c];

	for (sig_ot_it = &signals[sig_ot], ot = overtones; ot != ot_end; ++ot){
		*sig_ot_it++ = ot->ratio;
		*sig_ot_it++ = ot->phase;

		for (c = 0; c < frag->n_channels; ++c)
			*sig_ot_it++ = ot->levels.obj[c];
	}

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, sig_n, frag->rate))
		return -1;

	i = 0;

	while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig.len; ++i, ++j) {
			const double f = sig.vectors[sig_freq].data[j];
			const double max_ratio = half_rate / f;
			const double ph = sig.vectors[sig_phase].data[j];
			const double t = ph + origin + (double)i / frag->rate;
			const double m = k * f;
			const struct splat_vector *otv = &sig.vectors[sig_ot];

			for (ot = overtones; ot != ot_end; ++ot) {
				const double ratio = (otv++)->data[j];
				const double ot_ph = (otv++)->data[j];
				double s;

				if (ratio >= max_ratio) {
					otv += frag->n_channels;
					continue;
				}

				s = sin(m * ratio * (t + ot_ph));

				for (c = 0; c < frag->n_channels; ++c) {
					double x, y;

					x = sig.vectors[sig_amp + c].data[j];
					y = (otv++)->data[j];
					frag->data[c][i] += s * x * y;
				}
			}
		}
	}

	splat_signal_free(&sig);
	PyMem_Free(signals);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

/* -- Fast functions -- */

#if defined(SPLAT_NEON)
static float32x4_t splat_sine_neon(float32x4_t x)
{
	const float32x4_t *sine_tableq = (float32x4_t *)splat_sine_table;
	const float32x4_t pi = vdupq_n_f32(M_PI);
	float32x4_t polyq[4];
	float32x4_t a;
	float32x4_t b;
	float32x4_t p;
	float32x4_t y;
	uint32x4_t neg;
	uint32x4_t m;

	/* a = fmod(x, M_PI) */
	a = vmulq_f32(x, vdupq_n_f32(1.0 / M_PI));
	a = vcvtq_f32_u32(vcvtq_u32_f32(a));
	a = vmlsq_f32(x, a, pi);

	/* if a < 0 then a += M_PI */
	neg = vcltq_f32(a, vdupq_n_f32(0.0));
	neg = vandq_u32(neg, (uint32x4_t)pi);
	a = vaddq_f32(a, (float32x4_t)neg);

	/* m = int(a * table_len / M_PI) */
	m = vcvtq_u32_f32(vmulq_f32(a, splat_neon_sine_step));

	/* polyq1..4 = transpose(table[m]) */
	/* y = p0 + (x * p1) + (x^2 * p2) + (x^3 * p3) */
	polyq[0] = sine_tableq[m[0]];
	polyq[1] = sine_tableq[m[1]];
	*((float32x4x2_t *)&polyq[0]) = vtrnq_f32(polyq[0], polyq[1]);
	polyq[2] = sine_tableq[m[2]];
	polyq[3] = sine_tableq[m[3]];
	*((float32x4x2_t *)&polyq[2]) = vtrnq_f32(polyq[2], polyq[3]);
	y = vcombine_f32(vget_low_f32(polyq[0]), vget_low_f32(polyq[2]));
	p = vcombine_f32(vget_low_f32(polyq[1]), vget_low_f32(polyq[3]));
	y = vmlaq_f32(y, a, p);
	b = vmulq_f32(a, a);
	p = vcombine_f32(vget_high_f32(polyq[0]), vget_high_f32(polyq[2]));
	y = vmlaq_f32(y, b, p);
	b = vmulq_f32(b, a);
	p = vcombine_f32(vget_high_f32(polyq[1]), vget_high_f32(polyq[3]));
	y = vmlaq_f32(y, b, p);

	/* b = fmod(x, 2 * M_PI) */
	b = vmulq_f32(x, vdupq_n_f32(1.0 / (2.0 * M_PI)));
	b = vcvtq_f32_u32(vcvtq_u32_f32(b));
	b = vmlsq_f32(x, b, vdupq_n_f32(M_PI * 2.0));

	/* neg = 2 if b < M_PI */
	neg = vcgtq_f32(b, pi);
	neg = vshrq_n_u32(neg, 31);
	neg = vshlq_n_u32(neg, 1);

	/* b = -1 if M_PI < fmod(x, 2 * M_PI) < 2 * M_PI else 1 */
	b = vsubq_f32(vdupq_n_f32(1.0), vcvtq_f32_u32(neg));

	/* result = y * b */
	return vmulq_f32(y, b);
}
#elif defined(SPLAT_SSE)
static __m128 splat_sine_sse(__m128 x)
{
	__m128 xpi;
	__m128 a;
	__m128 b;
	__m128i m;
	uint32_t *mf = (uint32_t *)&m;
	__m128 poly2;
	__m128 poly3;
	__m128 poly4;
	__m128 neg;
	__m128 y;

	/* a = fmod(x, M_PI) */
	xpi = _mm_div_ps(x, splat_sse_pi);
	a = _mm_cvtepi32_ps(_mm_cvttps_epi32(xpi));
	a = _mm_sub_ps(x, _mm_mul_ps(a, splat_sse_pi));

	/* if a < 0 then a += M_PI (work around rounding issue) */
	neg = _mm_cmplt_ps(a, splat_sse_zero);
	neg = _mm_and_ps(neg, splat_sse_pi);
	a = _mm_add_ps(a, neg);

	/* if a > M_PI then a -= M_PI (work around rounding issue) */
	neg = _mm_cmpgt_ps(a, splat_sse_pi);
	neg = _mm_and_ps(neg, splat_sse_pi);
	a = _mm_sub_ps(a, neg);

	/* m = int(a * table_len / M_PI) */
	m = _mm_cvttps_epi32(_mm_mul_ps(a, splat_sse_sine_step));

	/* poly1..4 = transpose(table[m]) */
	y = _mm_load_ps(splat_sine_table[mf[0]].coef);
	poly2 = _mm_load_ps(splat_sine_table[mf[1]].coef);
	poly3 = _mm_load_ps(splat_sine_table[mf[2]].coef);
	poly4 = _mm_load_ps(splat_sine_table[mf[3]].coef);
	_MM_TRANSPOSE4_PS(y, poly2, poly3, poly4);

	/* y = p1 + (a * p2) + (a^2 * p3) + (a^3 * p4) */
	poly2 = _mm_mul_ps(poly2, a);
	y = _mm_add_ps(y, poly2);
	b = _mm_mul_ps(a, a);
	poly3 = _mm_mul_ps(poly3, b);
	y = _mm_add_ps(y, poly3);
	b = _mm_mul_ps(b, a);
	poly4 = _mm_mul_ps(poly4, b);
	y = _mm_add_ps(y, poly4);

	/* b = fmod(x, 2 * M_PI) */
	b = _mm_div_ps(x, splat_sse_pi2);
	b = _mm_cvtepi32_ps(_mm_cvttps_epi32(b));
	b = _mm_sub_ps(x, _mm_mul_ps(b, splat_sse_pi2));

	/* if b > M_PI then b = -1 else b = 1 */
	b = _mm_cmpgt_ps(b, splat_sse_pi);
	b = _mm_and_ps(b, splat_sse_two);
	b = _mm_sub_ps(splat_sse_one, b);

	/* result = y * b */
	return _mm_mul_ps(y, b);
}
#endif

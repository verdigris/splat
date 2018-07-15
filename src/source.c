/*
    Splat - source.c

    Copyright (C) 2015, 2017
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

#ifdef SPLAT_FAST
static sf_float_t splat_fast_sine(sf_float_t x);
#endif

/* -- sine source -- */

enum {
	SIG_SINE_FREQ = 0,
	SIG_SINE_PHASE,
	SIG_SINE_AMP,
} splat_sine_signal;

#ifdef SPLAT_FAST
void splat_sine_floats(struct splat_fragment *frag, const double *levels,
		       double freq, double phase)
{
	const sf_float_t phrate = sf_set(phase * frag->rate);
	const sf_float_t k = sf_set(2.0 * M_PI * freq / frag->rate);
	sf_float_t levelsq[SPLAT_MAX_CHANNELS];
	sf_float_t *out[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	for (c = 0; c < frag->n_channels; ++c) {
		levelsq[c] = sf_set(levels[c]);
		out[c] = (sf_float_t *)frag->channels[c].data;
	}

	for (i = 0; i < frag->length; i += 4) {
		sf_float_t x;
		sf_float_t y;

		/* x = 2 * M_PI * freq * (phrate + i) / frag->rate) */
		x = sf_set((float)i);
		x = sf_add(x, splat_fast_inc);
		x = sf_add(x, phrate);
		x = sf_mul(x, k);
		y = splat_fast_sine(x);

		for (c = 0; c < frag->n_channels; ++c)
			*out[c]++ = sf_mul(y, levelsq[c]);
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
			frag->channels[c].data[i] = s * levels[c];
	}
}
#endif

#ifdef SPLAT_FAST
static void _splat_sine_signals(struct splat_fragment *frag,
				struct splat_signal *sig, double origin)
{
	const sf_float_t k = sf_set(2.0 * M_PI);
#if defined(SPLAT_NEON)
	const sf_float_t rateq_inv = sf_set(1.0 / frag->rate);
#elif defined(SPLAT_SSE)
	const sf_float_t rateq = sf_set(frag->rate);
#endif
	const sf_float_t originq = sf_set(origin);
	sf_float_t *out[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i = 0;

	for (c = 0; c < frag->n_channels; ++c)
		out[c] = (sf_float_t *)frag->channels[c].data;

	while (splat_signal_next(sig) == SPLAT_SIGNAL_CONTINUE) {
		const sf_float_t *fq = (sf_float_t *)
			sig->vectors[SIG_SINE_FREQ].data;
		const sf_float_t *phq= (sf_float_t *)
			sig->vectors[SIG_SINE_PHASE].data;
		const sf_float_t *aq[SPLAT_MAX_CHANNELS];
		size_t j;

		for (c = 0; c < frag->n_channels; ++c)
			aq[c] = (sf_float_t *)
				sig->vectors[SIG_SINE_AMP + c].data;

		for (; c < SPLAT_MAX_CHANNELS; ++c)
			aq[c] = NULL;

		for (j = 0; j < sig->len; i += 4, j += 4) {
			sf_float_t x;
			sf_float_t f;
			sf_float_t y;

			x = sf_set((float)i);
			x = sf_add(x, splat_fast_inc);
#if defined(SPLAT_NEON)
			x = vmulq_f32(x, rateq_inv);
#elif defined(SPLAT_SSE)
			x = _mm_div_ps(x, rateq);
#endif
			x = sf_add(x, *phq++);
			x = sf_add(x, originq);
			f = sf_mul(*fq++, k);
			x = sf_mul(x, f);
			y = splat_fast_sine(x);

			for (c = 0; c < frag->n_channels; ++c)
				*out[c]++ = sf_mul(y, *aq[c]++);
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

				frag->channels[c].data[i] = s * a;
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
			frag->channels[c].data[i] = levels[c];
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

				frag->channels[c].data[i] = a * s;
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
			frag->channels[c].data[i] = (a[c] * t_rel) + b[c];
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

				frag->channels[c].data[i] = (a * t_rel) + b;
			}
		}
	}

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

/* -- overtones source -- */

enum {
	SIG_OT_FREQ = 0,
	SIG_OT_PHASE = 1,
	SIG_OT_AMP = 2,
};

#ifdef SPLAT_FAST
static void _splat_overtones_float(struct splat_fragment *frag,
				   double freq, double phase,
				   const struct splat_overtone *overtones,
				   const struct splat_overtone *ot_end)
{
	const sf_float_t k = sf_set(2 * M_PI * freq);
	const sf_float_t ph = sf_set(phase);
#if defined(SPLAT_NEON)
	const sf_float_t rate_inv = sf_set(1.0 / frag->rate);
#elif defined(SPLAT_SSE)
	const sf_float_t rateq = sf_set((float)frag->rate);
#endif
	const sf_float_t max_ratio = sf_set((frag->rate / 2.0) / freq);
	sf_float_t *outq[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t i;

	for (c = 0; c < frag->n_channels; ++c)
		outq[c] = (sf_float_t *)frag->channels[c].data;

	/* x = 2 * M_PI * freq * ratio * (ph + i/rate) */
	for (i = 0; i < frag->length; i += 4) {
		const struct splat_overtone *ot;
		sf_float_t x0;

		x0 = sf_set((float)i);
		x0 = sf_add(x0, splat_fast_inc);
#if defined(SPLAT_NEON)
		x0 = vmulq_f32(x0, rate_inv);
#elif defined(SPLAT_SSE)
		x0 = _mm_div_ps(x0, rateq);
#endif
		x0 = sf_add(x0, ph);

		for (ot = overtones; ot != ot_end; ++ot) {
			const sf_float_t *lvlq = ot->levels.flq;
			sf_float_t xr;
			sf_float_t yr;
			sf_mask_t clip;

			clip = sf_lt(ot->fl_ratioq, max_ratio);
			xr = sf_add(x0, ot->fl_phaseq);
			xr = sf_mul(xr, k);
			xr = sf_mul(xr, ot->fl_ratioq);
			yr = splat_fast_sine(xr);
			yr = sf_and(yr, clip);

			for (c = 0; c < frag->n_channels; ++c) {
				const sf_float_t y = sf_mul(yr, lvlq[c]);
				*outq[c] = sf_add(*outq[c], y);
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
				frag->channels[c].data[i] +=
					s * ot->levels.fl[c];
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

#ifdef SPLAT_FAST
		for (c = 0; c < frag->n_channels; ++c)
			ot->levels.flq[c] = sf_set(ot->levels.fl[c]);
#endif
	}

	_splat_overtones_float(frag, freq, phase, overtones, ot_end);
}

#ifdef SPLAT_FAST
static void _splat_overtones_mixed(struct splat_fragment *frag,
				   struct splat_signal *sig,
				   const struct splat_overtone *overtones,
				   size_t n, double origin)
{
	const sf_float_t k = sf_set(2.0 * M_PI);
#if defined(SPLAT_NEON)
	const sf_float_t rateq_inv = sf_set(1.0 / frag->rate);
#elif defined(SPLAT_SSE)
	const sf_float_t rateq = sf_set((float)frag->rate);
#endif
	const sf_float_t half_rate = sf_set(frag->rate / 2.0);
	const sf_float_t originq = sf_set(origin);
	const struct splat_overtone * const ot_end = &overtones[n];
	sf_float_t *out[SPLAT_MAX_CHANNELS];
	size_t i = 0;
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c)
		out[c] = (sf_float_t *)frag->channels[c].data;

	while (splat_signal_next(sig) == SPLAT_SIGNAL_CONTINUE) {
		const sf_float_t *fq = (sf_float_t *)
			sig->vectors[SIG_OT_FREQ].data;
		const sf_float_t *phq = (sf_float_t *)
			sig->vectors[SIG_OT_PHASE].data;
		const sf_float_t *amq[SPLAT_MAX_CHANNELS];
		size_t j;

		for (c = 0; c < frag->n_channels; ++c)
			amq[c] = (sf_float_t *)
				sig->vectors[SIG_OT_AMP + c].data;

		for (; c < SPLAT_MAX_CHANNELS; ++c)
			amq[c] = NULL;

		for (j = 0; j < sig->len; i += 4, j += 4) {
			const struct splat_overtone *ot;
			const sf_float_t f = *fq++;
#if defined(SPLAT_NEON)
			const sf_float_t max_ratio =
				vmulq_f32(half_rate, vrecpeq_f32(f));
#elif defined(SPLAT_SSE)
			const sf_float_t max_ratio = _mm_div_ps(half_rate, f);
#endif
			const sf_float_t fk = sf_mul(f, k);
			const sf_float_t pho = sf_add(*phq++, originq);
			sf_float_t aq[SPLAT_MAX_CHANNELS];
			sf_float_t y[SPLAT_MAX_CHANNELS];
			sf_float_t x0;

			for (c = 0; c < frag->n_channels; ++c) {
				aq[c] = *amq[c]++;
				y[c] = sf_zero();
			}

			for (; c < SPLAT_MAX_CHANNELS; ++c) {
				aq[c] = sf_zero();
				y[c] = sf_zero();
			}

			x0 = sf_set((float)i);
			x0 = sf_add(x0, splat_fast_inc);
#if defined(SPLAT_NEON)
			x0 = vmulq_f32(x0, rateq_inv);
#elif defined(SPLAT_SSE)
			x0 = _mm_div_ps(x0, rateq);
#endif
			x0 = sf_add(x0, pho);

			for (ot = overtones; ot != ot_end; ++ot) {
				const sf_float_t *lvlq = ot->levels.flq;
				sf_float_t x;
				sf_float_t yr;
				sf_mask_t clip;

				clip = sf_lt(ot->fl_ratioq, max_ratio);
				x = sf_add(x0, ot->fl_phaseq);
				x = sf_mul(x, ot->fl_ratioq);
				x = sf_mul(x, fk);
				yr = splat_fast_sine(x);

				for (c = 0; c < frag->n_channels; ++c) {
					sf_float_t y0;

					y0 = sf_mul(yr, aq[c]);
					y0 = sf_mul(y0, lvlq[c]);
					y0 = sf_and(y0, clip);
					y[c] = sf_add(y[c], y0);
				}
			}

			for (c = 0; c < frag->n_channels; ++c)
				*out[c]++ = y[c];
		}
	}
}
#else
static void _splat_overtones_mixed(struct splat_fragment *frag,
				   struct splat_signal *sig,
				   const struct splat_overtone *overtones,
				   size_t n, double origin)
{
	const double k = 2 * M_PI;
	const double half_rate = frag->rate / 2;
	const struct splat_overtone * const ot_end = &overtones[n];
	size_t i = 0;

	while (splat_signal_next(sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig->len; ++i, ++j) {
			const struct splat_overtone *ot;
			const double f = sig->vectors[SIG_OT_FREQ].data[j];
			const double max_ratio = half_rate / f;
			const double m = k * f;
			const double ph = sig->vectors[SIG_OT_PHASE].data[j];
			const double t = ph + origin + (double)i / frag->rate;

			for (ot = overtones; ot != ot_end; ++ot) {
				unsigned c;
				double s;

				if (ot->fl_ratio >= max_ratio)
					continue;

				s = sin(m * ot->fl_ratio * (t + ot->fl_phase));

				for (c = 0; c < frag->n_channels; ++c) {
					double x;

					x = sig->vectors[SIG_OT_AMP+c].data[j];
					x *= ot->levels.fl[c];
					frag->channels[c].data[i] += s * x;
				}
			}
		}
	}
}
#endif

int splat_overtones_mixed(struct splat_fragment *frag, PyObject **levels,
			  PyObject *freq, PyObject *phase,
			  struct splat_overtone *overtones, Py_ssize_t n,
			  double origin)
{
	PyObject *signals[SIG_OT_AMP + SPLAT_MAX_CHANNELS];
	struct splat_signal sig;
	unsigned c;

	signals[SIG_OT_FREQ] = freq;
	signals[SIG_OT_PHASE] = phase;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_OT_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, (SIG_OT_AMP + frag->n_channels),
			      frag->rate))
		return -1;

	_splat_overtones_mixed(frag, &sig, overtones, n, origin);

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

#ifdef SPLAT_FAST
static void _splat_overtones_signal(struct splat_fragment *frag,
				    struct splat_signal *sig,
				    const struct splat_overtone *overtones,
				    size_t n, double origin)
{
	const sf_float_t k = sf_set(2.0 * M_PI);
#if defined(SPLAT_NEON)
	const sf_float_t rateq_inv = sf_set(1.0 / frag->rate);
#elif defined(SPLAT_SSE)
	const sf_float_t rateq = sf_set((float)frag->rate);
#endif
	const sf_float_t half_rate = sf_set(frag->rate / 2.0);
	const sf_float_t originq = sf_set(origin);
	const struct splat_overtone * const ot_end = &overtones[n];
	const size_t sig_ot = SIG_OT_AMP + frag->n_channels;
	sf_float_t *out[SPLAT_MAX_CHANNELS];
	size_t i = 0;
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c)
		out[c] = (sf_float_t *)frag->channels[c].data;

	while (splat_signal_next(sig) == SPLAT_SIGNAL_CONTINUE) {
		const sf_float_t *fq = (sf_float_t *)
			sig->vectors[SIG_OT_FREQ].data;
		const sf_float_t *phq = (sf_float_t *)
			sig->vectors[SIG_OT_PHASE].data;
		const sf_float_t *amq[SPLAT_MAX_CHANNELS];
		size_t j;

		for (c = 0; c < frag->n_channels; ++c)
			amq[c] = (sf_float_t *)
				sig->vectors[SIG_OT_AMP + c].data;

		for (; c < SPLAT_MAX_CHANNELS; ++c)
			amq[c] = NULL;

		for (j = 0; j < sig->len; i += 4, j += 4) {
			const struct splat_vector *otv = &sig->vectors[sig_ot];
			const struct splat_overtone *ot;
			const sf_float_t f = *fq++;
#if defined(SPLAT_NEON)
			const sf_float_t max_ratio =
				vmulq_f32(half_rate, vrecpeq_f32(f));
#elif defined(SPLAT_SSE)
			const sf_float_t max_ratio =
				_mm_div_ps(half_rate, f);
#endif
			const sf_float_t fk = sf_mul(f, k);
			const sf_float_t pho = sf_add(*phq++, originq);
			sf_float_t aq[SPLAT_MAX_CHANNELS];
			sf_float_t y[SPLAT_MAX_CHANNELS];
			sf_float_t x0;

			for (c = 0; c < frag->n_channels; ++c) {
				aq[c] = *amq[c]++;
				y[c] = sf_zero();
			}

			for (; c < SPLAT_MAX_CHANNELS; ++c) {
				aq[c] = sf_zero();
				y[c] = sf_zero();
			}

			x0 = sf_set((float)i);
			x0 = sf_add(x0, splat_fast_inc);
#if defined(SPLAT_NEON)
			x0 = vmulq_f32(x0, rateq_inv);
#elif defined(SPLAT_SSE)
			x0 = _mm_div_ps(x0, rateq);
#endif
			x0 = sf_add(x0, pho);

			for (ot = overtones; ot != ot_end; ++ot) {
				const sf_float_t ot_r =
					*((sf_float_t *)&(otv++)->data[j]);
				const sf_float_t ot_ph =
					*((sf_float_t *)&(otv++)->data[j]);
				sf_mask_t clip;
				sf_float_t x;
				sf_float_t yr;

				clip = sf_lt(ot_r, max_ratio);
				x = sf_add(x0, ot_ph);
				x = sf_mul(x, ot_r);
				x = sf_mul(x, fk);
				yr = splat_fast_sine(x);

				for (c = 0; c < frag->n_channels; ++c) {
					sf_float_t g;
					sf_float_t y0;

					g = *((sf_float_t*)&(otv++)->data[j]);
					y0 = sf_mul(yr, aq[0]);
					y0 = sf_mul(y0, g);
					y0 = sf_and(y0, clip);
					y[c] = sf_add(y[c], y0);
				}
			}

			for (c = 0; c < frag->n_channels; ++c)
				*out[c]++ = y[c];
		}
	}
}
#else
static void _splat_overtones_signal(struct splat_fragment *frag,
				    struct splat_signal *sig,
				    const struct splat_overtone *overtones,
				    size_t n, double origin)
{
	const double k = 2 * M_PI;
	const double half_rate = frag->rate / 2;
	const struct splat_overtone * const ot_end = &overtones[n];
	const size_t sig_ot = SIG_OT_AMP + frag->n_channels;
	size_t i = 0;

	while (splat_signal_next(sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig->len; ++i, ++j) {
			const double f = sig->vectors[SIG_OT_FREQ].data[j];
			const double max_ratio = half_rate / f;
			const double ph = sig->vectors[SIG_OT_PHASE].data[j];
			const double t = ph + origin + (double)i / frag->rate;
			const double m = k * f;
			const struct splat_vector *otv = &sig->vectors[sig_ot];
			const struct splat_overtone *ot;

			for (ot = overtones; ot != ot_end; ++ot) {
				const double ratio = (otv++)->data[j];
				const double ot_ph = (otv++)->data[j];
				unsigned c;
				double s;

				if (ratio >= max_ratio) {
					otv += frag->n_channels;
					continue;
				}

				s = sin(m * ratio * (t + ot_ph));

				for (c = 0; c < frag->n_channels; ++c) {
					double g;

					g = sig->vectors[SIG_OT_AMP+c].data[j];
					g *= (otv++)->data[j];
					frag->channels[c].data[i] += s * g;
				}
			}
		}
	}
}
#endif

int splat_overtones_signal(struct splat_fragment *frag, PyObject **levels,
			   PyObject *freq, PyObject *phase,
			   struct splat_overtone *overtones, Py_ssize_t n,
			   double origin)
{
	PyObject **signals;
	struct splat_signal sig;
	const struct splat_overtone *ot;
	const size_t sig_ot = SIG_OT_AMP + frag->n_channels;
	const struct splat_overtone * const ot_end = &overtones[n];
	PyObject **sig_ot_it;
	/* for each overtone: ratio, phase and levels */
	const size_t sig_n = sig_ot + (n * (2 + frag->n_channels));
	unsigned c;

	signals = PyMem_Malloc(sig_n * sizeof(PyObject *));

	if (signals == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	signals[SIG_OT_FREQ] = freq;
	signals[SIG_OT_PHASE] = phase;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_OT_AMP + c] = levels[c];

	for (sig_ot_it = &signals[sig_ot], ot = overtones; ot != ot_end; ++ot){
		*sig_ot_it++ = ot->ratio;
		*sig_ot_it++ = ot->phase;

		for (c = 0; c < frag->n_channels; ++c)
			*sig_ot_it++ = ot->levels.obj[c];
	}

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, sig_n, frag->rate))
		return -1;

	_splat_overtones_signal(frag, &sig, overtones, n, origin);

	splat_signal_free(&sig);
	PyMem_Free(signals);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

/* -- Fast functions -- */

#if defined(SPLAT_NEON)
static float32x4_t splat_fast_sine(float32x4_t x)
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

	/* m = int(x * table_len / M_PI) & table_mask */
	m = vcvtq_u32_f32(vmulq_f32(x, splat_fast_sine_step));
	m = vandq_u32(m, splat_fast_sine_mask);

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
static __m128 splat_fast_sine(__m128 x)
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
	m = _mm_cvttps_epi32(_mm_mul_ps(a, splat_fast_sine_step));

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

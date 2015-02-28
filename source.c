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

/* -- sine source -- */

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

int splat_sine_signals(struct splat_fragment *frag, PyObject **levels,
		       PyObject *freq, PyObject *phase, double origin)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_AMP,
	};
	static const double k = 2 * M_PI;
	struct splat_signal sig;
	PyObject *signals[SIG_AMP + SPLAT_MAX_CHANNELS];
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
			const double t = (double)i / frag->rate;
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double s = sin(k * f * (t + ph + origin));

			for (c = 0; c < frag->n_channels; ++c) {
				const double a =
					sig.vectors[SIG_AMP + c].data[j];

				frag->data[c][i] = s * dB2lin(a);
			}
		}
	}

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}


/* -- square source -- */

void splat_square_floats(struct splat_fragment *frag, const double *fl_pos,
			 double freq, double phase, double ratio)
{
	const double k = freq / frag->rate;
	const double ph0 = freq * phase;
	double fl_neg[SPLAT_MAX_CHANNELS];
	unsigned c;
	Py_ssize_t i;

	ratio = min(ratio, 1.0);
	ratio = max(ratio, 0.0);

	for (c = 0; c < frag->n_channels; ++c)
		fl_neg[c] = -fl_pos[c];

	for (i = 0; i < frag->length; ++i) {
		double n_periods;
		const double t_rel = modf(((i * k) + ph0), &n_periods);
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
			const double t = (double)i / frag->rate;
			const double tph = t + ph + origin;
			double n_periods;
			const double t_rel = modf((f * (tph)), &n_periods);
			double ratio = sig.vectors[SIG_RATIO].data[j];
			double s;

			ratio = min(ratio, 1.0);
			ratio = max(ratio, 0.0);
			s = (t_rel < ratio) ? 1.0 : -1.0;

			for (c = 0; c < frag->n_channels; ++c) {
				const double a =
					sig.vectors[SIG_AMP + c].data[j];

				frag->data[c][i] = dB2lin(a) * s;
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
	const double ph0 = freq * phase;
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
		const double t_rel = modf(((i * k) + ph0), &n_periods);

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
			const double t = (double)i / frag->rate;
			const double tph = t + ph + origin;
			double n_periods;
			const double t_rel = modf((f * (tph)), &n_periods);
			double ratio = sig.vectors[SIG_RATIO].data[j];

			ratio = min(ratio, 1.0);
			ratio = max(ratio, 0.0);

			for (c = 0; c < frag->n_channels; ++c) {
				const double l_log =
					sig.vectors[SIG_AMP + c].data[j];
				const double l = dB2lin(l_log);
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

void splat_overtones_float(struct splat_fragment *frag, const double *levels,
			   double freq, double phase,
			   struct splat_overtone *overtones, Py_ssize_t n)
{
	const double k = 2 * M_PI * freq;
	const double max_ratio = (frag->rate / freq) / 2;
	struct splat_overtone *ot;
	const struct splat_overtone *ot_end = &overtones[n];
	unsigned c;
	size_t i;

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

	for (i = 0; i < frag->length; ++i) {
		const double t = phase + ((double)i / frag->rate);

		for (ot = overtones; ot != ot_end; ++ot) {
			const double s =
				sin(k * ot->fl_ratio * (t + ot->fl_phase));

			for (c = 0; c < frag->n_channels; ++c)
				frag->data[c][i] += s * ot->levels.fl[c];
		}
	}
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
			const double t = (double)i / frag->rate;
			const double tph = t + ph + origin;

			for (ot = overtones; ot != ot_end; ++ot) {
				double s;

				if (ot->fl_ratio >= max_ratio)
					continue;

				s = sin(m * ot->fl_ratio *
					(tph + ot->fl_phase));

				for (c = 0; c < frag->n_channels; ++c) {
					double x;

					x = sig.vectors[SIG_AMP + c].data[j];
					x = dB2lin(x);
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

	{
		PyObject **sig_ot_it = &signals[sig_ot];

		for (ot = overtones; ot != ot_end; ++ot) {
			*sig_ot_it++ = ot->ratio;
			*sig_ot_it++ = ot->phase;

			for (c = 0; c < frag->n_channels; ++c)
				*sig_ot_it++ = ot->levels.obj[c];
		}
	}

	if (splat_signal_init(&sig, frag->length, (origin * frag->rate),
			      signals, sig_n, frag->rate))
		return -1;

	i = 0;

	while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig.len; ++i, ++j) {
			const double t = (double)i / frag->rate;
			const double f = sig.vectors[sig_freq].data[j];
			const double max_ratio = half_rate / f;
			const double ph = sig.vectors[sig_phase].data[j];
			const double m = k * f;
			const double tph = t + ph + origin;
			const struct splat_vector *otv = &sig.vectors[sig_ot];

			for (ot = overtones; ot != ot_end; ++ot) {
				const double ratio = (otv++)->data[j];
				const double ot_ph = (otv++)->data[j];
				double s;

				if (ratio >= max_ratio) {
					otv += frag->n_channels;
					continue;
				}

				s = sin(m * ratio * (tph + ot_ph));

				for (c = 0; c < frag->n_channels; ++c) {
					double x, y;

					x = sig.vectors[sig_amp + c].data[j];
					y = (otv++)->data[j];
					frag->data[c][i] +=
						s * dB2lin(x) * dB2lin(y);
				}
			}
		}
	}

	splat_signal_free(&sig);
	PyMem_Free(signals);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}


/*
    Splat - frag.c

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

static void *splat_frag_realloc(void *data, size_t cur_size, size_t new_size)
{
	void *new_data;

#ifdef SPLAT_FAST
	if (posix_memalign(&new_data, 16, new_size))
		return NULL;

	if (data != NULL) {
		memcpy(new_data, data, min(cur_size, new_size));
		free(data);
	}
#else
	if (data == NULL)
		new_data = PyMem_Malloc(new_size);
	else
		new_data = PyMem_Realloc(data, new_size);
#endif

	return new_data;
}

int splat_frag_init(struct splat_fragment *frag, unsigned n_channels,
		    unsigned rate, size_t length, const char *name)
{
	size_t data_size;
	unsigned c;

#ifdef SPLAT_FAST
	length = splat_round4(length);
#endif
	data_size = length * sizeof(sample_t);

	for (c = 0; c < n_channels; ++c) {
		if (!length) {
			frag->data[c] = NULL;
		} else {
#ifdef SPLAT_FAST
			void *data;

			if (posix_memalign(&data, 16, data_size))
				frag->data[c] = NULL;
			else
				frag->data[c] = data;
#else
			frag->data[c] = PyMem_Malloc(data_size);
#endif
			if (frag->data[c] == NULL) {
				PyErr_NoMemory();
				return -1;
			}

			memset(frag->data[c], 0, data_size);
		}
	}

	if (name == NULL)
		frag->name = NULL;
	else if (splat_frag_set_name(frag, name))
		return -1;

	frag->n_channels = n_channels;
	frag->rate = rate;
	frag->length = length;

	return 0;
}

void splat_frag_free(struct splat_fragment *frag)
{
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c)
#ifdef SPLAT_FAST
		free(frag->data[c]);
#else
		PyMem_Free(frag->data[c]);
#endif

	if (frag->name != NULL)
		free(frag->name);
}

int splat_frag_set_name(struct splat_fragment *frag, const char *name)
{
	if (frag->name != NULL)
		free(frag->name);

	frag->name = strdup(name);

	if (frag->name == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	return 0;
}

int splat_frag_resize(struct splat_fragment *frag, size_t length)
{
	size_t start = frag->length * sizeof(sample_t);
	size_t size;
	ssize_t extra;
	unsigned c;

	if (length == frag->length)
		return 0;

	start = frag->length * sizeof(sample_t);
#if SPLAT_FAST
	length = splat_round4(length);
#endif
	size = length * sizeof(sample_t);
	extra = size - start;

	for (c = 0; c < frag->n_channels; ++c) {
		frag->data[c] = splat_frag_realloc(frag->data[c], start, size);

		if (frag->data[c] == NULL) {
			PyErr_NoMemory();
			return -1;
		}

		if (extra > 0)
			memset(&frag->data[c][frag->length], 0, extra);
	}

	frag->length = length;

	return 0;
}

#ifdef SPLAT_FAST
static void splat_frag_mix_floats(struct splat_fragment *frag,
				  const struct splat_fragment *in,
				  size_t offset, size_t start, size_t length,
				  const double *levels, int zero_dB)
{
	unsigned c;

	start = splat_mask4(start);
	offset = splat_mask4(offset);
	length /= 4;

	for (c = 0; c < frag->n_channels; ++c) {
		const sf_float_t *src = (sf_float_t *)&in->data[c][start];
		sf_float_t *dst = (sf_float_t *)&frag->data[c][offset];
		size_t i = length;

		if (zero_dB) {
			while (i--) {
				*dst = sf_add(*dst, *src++);
				dst++;
			}
		} else {
			const sf_float_t gain = sf_set(levels[c]);

			while (i--) {
#if defined(SPLAT_NEON)
				*dst = vmlaq_f32(*dst, *src++, gain);
#else
				sf_float_t s = *src++;

				s = sf_mul(s, gain);
				*dst = sf_add(*dst, s);
#endif
				dst++;
			}
		}
	}
}
#else
static void splat_frag_mix_floats(struct splat_fragment *frag,
				  const struct splat_fragment *incoming,
				  size_t offset, size_t start, size_t length,
				  const double *levels, int zero_dB)
{
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c) {
		const sample_t *src = &incoming->data[c][start];
		sample_t *dst =  &frag->data[c][offset];
		size_t i = length;

		if (zero_dB) {
			while (i--)
				*dst++ += *src++;
		} else {
			const double g = levels[c];

			while (i--)
				*dst++ += g * (*src++);
		}
	}
}
#endif

static int splat_frag_mix_signals(struct splat_fragment *frag,
				  const struct splat_fragment *incoming,
				  size_t offset, size_t start, size_t length,
				  const struct splat_levels *levels)
{
	struct splat_signal sig;
	PyObject *signals[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t in;
	size_t i;

	for (c = 0; c < incoming->n_channels; ++c)
		signals[c] = levels->obj[c];

	if (splat_signal_init(&sig, length, offset, signals, frag->n_channels,
			      frag->rate))
		return -1;

	in = sig.cur;
	i = 0;

	while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig.len; ++i, ++j, ++in) {
			for (c = 0; c < incoming->n_channels; ++c) {
				double a = sig.vectors[c].data[j];

				frag->data[c][i] += incoming->data[c][in] * a;
			}
		}
	}

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

int splat_frag_mix(struct splat_fragment *frag,
		   const struct splat_fragment *incoming,
		   const struct splat_levels *levels,
		   size_t length, double offset, double skip, int zero_dB)
{
	size_t offset_sample;
	size_t skip_sample;
	size_t total_length;

	offset_sample = offset * frag->rate;
	offset_sample = max(offset_sample, 0);
	skip_sample = skip * frag->rate;
	skip_sample = minmax(skip_sample, 0, incoming->length);
	length = minmax(length, 0, (incoming->length - skip_sample));
	total_length = offset_sample + length;

	if (splat_frag_grow(frag, total_length))
		return -1;

	if (levels->all_floats)
		splat_frag_mix_floats(frag, incoming, offset_sample,
				      skip_sample, length, levels->fl,
				      zero_dB);
	else if (splat_frag_mix_signals(frag, incoming, offset_sample,
					skip_sample, length, levels))
		return -1;

	return 0;
}

int splat_frag_sample_number(size_t *val, long min_val, long max_val,
			     PyObject *obj)
{
	long tmp_val;

	if (!PyInt_Check(obj)) {
		PyErr_SetString(PyExc_TypeError,
				"sample number must be an integer");
		return -1;
	}

	tmp_val = min(PyInt_AS_LONG(obj), max_val);
	*val = max(tmp_val, min_val);

	return 0;
}

void splat_frag_get_peak(const struct splat_fragment *frag,
			 struct splat_peak *chan_peak,
			 struct splat_peak *frag_peak, int do_avg)
{
	unsigned c;

	frag_peak->avg = 0.0;
	frag_peak->max = -1.0;
	frag_peak->min = 1.0;
	frag_peak->peak = 0.0;

	for (c = 0; c < frag->n_channels; ++c) {
		sample_t * const chan_data = frag->data[c];
		const sample_t * const end = &chan_data[frag->length];
		const sample_t *it;
		double avg = 0.0;
		double max = -1.0;
		double min = 1.0;

		for (it = chan_data; it != end; ++it) {
			if (do_avg)
				avg += *it / frag->length;

			if (*it > max)
				max = *it;
			else if (*it < min)
				min = *it;
		}

		chan_peak[c].avg = avg;
		chan_peak[c].max = max;
		chan_peak[c].min = min;

		if (do_avg)
			frag_peak->avg += avg / frag->n_channels;

		if (max > frag_peak->max)
			frag_peak->max = max;

		if (min < frag_peak->min)
			frag_peak->min = min;

		min = fabsf(min);
		max = fabsf(max);
		chan_peak[c].peak = (min > max) ? min : max;

		if (frag_peak->peak < chan_peak[c].peak)
			frag_peak->peak = chan_peak[c].peak;
	}
}

void splat_frag_normalize(struct splat_fragment *frag, double level_dB,
			  int do_zero)
{
	const double level = dB2lin(level_dB);
	struct splat_peak chan_peak[SPLAT_MAX_CHANNELS];
	struct splat_peak frag_peak;
	unsigned c;
	double gain;

	splat_frag_get_peak(frag, chan_peak, &frag_peak, do_zero);

	if (do_zero) {
		double offset = 0.0;

		for (c = 0; c < frag->n_channels; ++c) {
			const double avg = chan_peak[c].avg;

			if (fabs(offset) < fabs(avg))
				offset = avg;
		}

		gain = level / (frag_peak.peak + fabs(offset));
	} else {
		gain = level / frag_peak.peak;
	}

	if ((1.0 < gain) && (gain < 1.001)) {
		int zero;

		if (!do_zero)
			return;

		for (c = 0, zero = 1; c < frag->n_channels && zero; ++c)
			if (fabs(chan_peak[c].avg) > 0.001)
				zero = 0;

		if (zero)
			return;
	}

	for (c = 0; c < frag->n_channels; ++c) {
		const double chan_avg = chan_peak[c].avg;
		sample_t * const chan_data = frag->data[c];
		const sample_t * const end = &chan_data[frag->length];
		sample_t *it;

		for (it = chan_data; it != end; ++it) {
			*it -= chan_avg;
			*it *= gain;
		}
	}
}

static void splat_frag_amp_floats(struct splat_fragment *frag,
				  const double *gains)
{
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c) {
		const double g = gains[c];
		size_t i;

		if (g == 1.0)
			continue;

		for (i = 0; i < frag->length; ++i)
			frag->data[c][i] *= g;
	}
}

static int splat_frag_amp_signals(struct splat_fragment *frag,
				  const struct splat_levels *gains)
{
	struct splat_signal sig;
	PyObject *signals[SPLAT_MAX_CHANNELS];
	unsigned c;
	size_t in;
	size_t i;

	for (c = 0; c < frag->n_channels; ++c)
		signals[c] = gains->obj[c];

	if (splat_signal_init(&sig, frag->length, 0.0, signals,
			      frag->n_channels, frag->rate))
		return -1;

	in = sig.cur;
	i = 0;

	while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
		size_t j;

		for (j = 0; j < sig.len; ++i, ++j, ++in)
			for (c = 0; c < frag->n_channels; ++c)
				frag->data[c][i] *= sig.vectors[c].data[j];
	}

	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;
}

int splat_frag_amp(struct splat_fragment *frag, struct splat_levels *gains)
{
	if (gains->all_floats)
		splat_frag_amp_floats(frag, gains->fl);
	else if (splat_frag_amp_signals(frag, gains))
		return -1;

	return 0;
}

void splat_frag_lin2dB(struct splat_fragment *frag)
{
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c) {
		size_t i;

		for (i = 0; i < frag->length; ++i)
			frag->data[c][i] = lin2dB(frag->data[c][i]);
	}
}

void splat_frag_dB2lin(struct splat_fragment *frag)
{
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c) {
		size_t i;

		for (i = 0; i < frag->length; ++i)
			frag->data[c][i] = dB2lin(frag->data[c][i]);
	}
}

int splat_frag_offset(struct splat_fragment *frag, PyObject *offset_obj,
		      double start)
{
	size_t i;
	unsigned c;

	if (PyFloat_Check(offset_obj)) {
		const double offset_float = PyFloat_AS_DOUBLE(offset_obj);

		for (c = 0; c < frag->n_channels; ++c) {
			for (i = 0; i < frag->length; ++i)
				frag->data[c][i] += offset_float;
		}
	} else {
		struct splat_signal sig;

		if (splat_signal_init(&sig, frag->length, (start * frag->rate),
				      &offset_obj, 1, frag->rate))
			return -1;

		i = 0;

		while (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE) {
			size_t j;

			for (j = 0; j < sig.len; ++i, ++j) {
				const double value = sig.vectors[0].data[j];

				for (c = 0; c < frag->n_channels; ++c)
					frag->data[c][i] += value;
			}
		}

		splat_signal_free(&sig);

		if (sig.stat == SPLAT_SIGNAL_ERROR)
			return -1;
	}

	return 0;
}

/*
  This resampling method uses quadratic interpolation.  The k coefficients are:

  y = k0 + k1 * x + k2 * x * x

  (1.1) y0 = k0 + k1 * x0 + k2 * x0 * x0
  (1.2) y1 = k0 + k1 * x1 + k2 * x1 * x1
  (1.3) y2 = k0 + k1 * x2 + k2 * x2 * x2

  x0 = x1 - 1
  x2 = x1 + 1

  (2.1) y0 = k0 + k1 * (x1 - 1) + k2 * (x1 * x1 - 2 * x1 + 1)
  (2.2) y1 = k0 + k1 *  x1      + k2 *  x1 * x1
  (2.3) y2 = k0 + k1 * (x1 + 1) + k2 * (x1 * x1 + 2 * x1 + 1)

  (4.2): (2.2)
    k0 = y1 - k1 * x1 - k2 * x1 * x1

  (4.1): (2.1), (4.2)
    y0 = y1 - k1 * x1 - k2 * x1 * x1 +
         k1 * (x1 - 1) + k2 * (x1 * x1 - 2 * x1 + 1)
       = y1 + k1 * (-x1 + x1 - 1) + k2 * (-x1 * x1 + x1 * x1 - 2 * x1 + 1)
       = y1 + k1 * (-1) + k2 * (-2 * x1 + 1)
       = y1 - k1 + k2 * (1 - 2 * x1)
   -y0 = -y1 + k1 - k2 * (1 - 2 * x1)
    k1 = y1 - y0 + k2 * (1 - 2 * x1)
    --------------------------------

  (5.1): (4.1), (4.2)
    k0 = y1 - (y1 - y0 + k2 * (1 - 2 * x1)) * x1 - k2 * x1 * x1
       = y1 - x1 * (y1 - y0) - k2 * (1 - 2 * x1) * x1 - k2 * x1 * x1
       = y1 - x1 * (y1 - y0) + k2 * (-x1 + 2 * x1 * x1 - x1 * x1)
    k0 = y1 + x1 * (y0 - y1) + k2 * (x1 * x1 - x1)
    ----------------------------------------------

  (5.3): (5.1), (4.2), (2.3)
    y2 = y1 - k1 * x1 - k2 * x1 * x1 +
         (y1 - y0 + k2 * (1 - 2 * x1)) * (x1 + 1) +
         k2 * (x1 * x1 + 2 * x1 + 1)
       = y1 - (y1 - y0 + k2 * (1 - 2 * x1)) * x1 - k2 * x1 * x1 +
         (y1 - y0 + k2 * (1 - 2 * x1)) * (x1 + 1) +
         k2 * (x1 * x1 + 2 * x1 + 1)
       = y1 + x1 * (y0 - y1) + k2 * (x1 * (2 * x1 - 1) - x1 * x1) +
         (x1 + 1) * (y1 - y0) + k2 * (1 - 2 * x1) * (x1 + 1) +
         k2 * (x1 * x1 + 2 * x1 + 1)
       = y1 + x1 * (y0 - y1) + (x1 + 1) * (y1 - y0) +
         k2 * (x1 * (2 * x1 - 1) - x1 * x1 +
         (1 - 2 * x1) * (x1 + 1) +
         (x1 * x1 + 2 * x1 + 1))
       = y1 * 2 - y0 +
         k2 * (x1 * x1 - x1 +
         1 - 2 * x1 + x1 - 2 * x1 * x1 +
         x1 * x1 + 2 * x1 + 1)
       = y1 * 2 - y0 +
         k2 * (2 + x1 * (-1 + 1 + 2 - 2) + (x1 * x1) * (1 - 2 + 1))
       = y1 * 2 - y0 + k2 * 2
    k2 = (y2 + y0 - 2 * y1) / 2
    ---------------------------
*/

int splat_frag_resample(struct splat_fragment *frag,
			const struct splat_fragment *old_frag,
			unsigned rate, double time_ratio)
{
	const double ratio = time_ratio * (double)rate / frag->rate;
	const size_t max_j = old_frag->length - 2;
	unsigned c;

	if (time_ratio <= 0.0) {
		PyErr_SetString(PyExc_ValueError,
				"resample time ratio must be positive");
		return -1;
	}

	if (splat_frag_resize(frag, frag->length * ratio))
		return -1;

	frag->rate = rate;

	if (!frag->length || !ratio)
		return 0;

	for (c = 0; c < frag->n_channels; ++c) {
		sample_t *to = frag->data[c];
		size_t i;

		*to++ = old_frag->data[c][0];

		for (i = 1; i < frag->length; ++i) {
			const double x = i / ratio;
			size_t j = x + 0.5;
			double x1;
			double y0, y1, y2;
			double k0, k1, k2;

			j = minmax(j, 1, max_j);
			x1 = j;
			y0 = old_frag->data[c][j - 1];
			y1 = old_frag->data[c][j++];
			y2 = old_frag->data[c][j];
			k2 = (y2 + y0 - 2 * y1) / 2;
			k1 = y1 - y0 + k2 * (1 - 2 * x1);
			k0 = y1 + x1 * (y0 - y1) + k2 * (x1 * x1 - x1);
			*to++ = k0 + x * k1 + x * x * k2;
		}
	}

	return 0;
}

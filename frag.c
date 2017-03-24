/*
    Splat - frag.c

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

static void *splat_channel_realloc(void *data, size_t cur_size,
				   size_t new_size)
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

static int splat_channel_resize_mem(struct splat_channel *chan, size_t length)
{
	const size_t start = chan->length * sizeof(sample_t);
	const size_t size = length * sizeof(sample_t);
	const ssize_t extra = size - start;

	chan->data = splat_channel_realloc(chan->data, start, size);

	if (chan->data == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	if (extra > 0)
		memset(&chan->data[chan->length], 0, extra);

	chan->length = length;

	return 0;
}

static void splat_channel_free_mem(struct splat_channel *chan)
{
#ifdef SPLAT_FAST
	free(chan->data);
#else
	PyMem_Free(chan->data);
#endif
}

static int splat_channel_init_mem(struct splat_channel *chan, size_t length)
{
	size_t data_size;
#ifdef SPLAT_FAST
	void *data;
#endif

	chan->length = length;
	chan->free = splat_channel_free_mem;
	chan->resize = splat_channel_resize_mem;
	chan->mmap.ptr = NULL;

	if (!length) {
		chan->data = NULL;
		return 0;
	}

	data_size = length * sizeof(sample_t);

#ifdef SPLAT_FAST
	if (posix_memalign(&data, 16, data_size))
		chan->data = NULL;
	else
		chan->data = data;
#else
	chan->data = PyMem_Malloc(data_size);
#endif

	if (chan->data == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	memset(chan->data, 0, data_size);

	return 0;
}

static int splat_channel_resize_mmap(struct splat_channel *chan, size_t length)
{
	const size_t size = length * sizeof(sample_t);

	if (splat_mmap_remap(&chan->mmap, size)) {
		PyErr_SetString(PyExc_SystemError, "mmap failed");
		return -1;
	}

	chan->data = chan->mmap.ptr;
	chan->length = length;

	return 0;
}

static void splat_channel_free_mmap(struct splat_channel *chan)
{
	splat_mmap_free(&chan->mmap);
}

static int splat_channel_init_mmap(struct splat_channel *chan, size_t length,
				   const char *mmap_path, int temp_mmap)
{
	chan->free = splat_channel_free_mmap;
	chan->resize = splat_channel_resize_mmap;

	if (splat_mmap_init(&chan->mmap, mmap_path)) {
		PyErr_SetString(PyExc_SystemError, "mmap init error");
		return -1;
	}

	if (chan->resize(chan, length))
		return -1;

	chan->mmap.persist = temp_mmap ? 0 : 1;
	chan->data = chan->mmap.ptr;

	return 0;
}

static int splat_frag_init_common(struct splat_fragment *frag,
				  unsigned n_channels, unsigned rate,
				  size_t length, const char *name)
{
	if (name == NULL)
		frag->name = NULL;
	else if (splat_frag_set_name(frag, name))
		return -1;

	frag->n_channels = n_channels;
	frag->rate = rate;
	frag->length = length;

	return 0;
}

int splat_frag_init(struct splat_fragment *frag, unsigned n_channels,
		    unsigned rate, size_t length, const char *name)
{
	unsigned c;

#if SPLAT_FAST
	length = splat_round4(length);
#endif

	for (c = 0; c < n_channels; ++c) {

		if (splat_channel_init_mem(&frag->channels[c], length))
			return -1;
	}

	frag->uses_mmap = 0;

	return splat_frag_init_common(frag, n_channels, rate, length, name);
}

int splat_frag_init_mmap(struct splat_fragment *frag, unsigned n_channels,
			 unsigned rate, size_t length, const char *name,
			 const char *new_path, const char **open_paths)
{
	char tmp_path[L_tmpnam];
	char *mmap_chan_path;
	size_t mmap_chan_path_len;
	unsigned c;

	frag->uses_mmap = 1;
	frag->temp_mmap = 0;

	if ((open_paths == NULL) || (open_paths[0] == NULL)) {
		if (new_path == NULL) {
			/* ToDo: use mkstemp instead */
			if (tmpnam_r(tmp_path) == NULL) {
				PyErr_SetString(PyExc_SystemError,
						"failed to create temp file");
				return -1;
			}
			new_path = tmp_path;
			frag->temp_mmap = 1;
		}

		mmap_chan_path_len = strlen(new_path) + 8;
		mmap_chan_path = PyMem_Malloc(mmap_chan_path_len);

		if (mmap_chan_path == NULL) {
			PyErr_NoMemory();
			return -1;
		}
	} else {
		mmap_chan_path = NULL;
	}

	for (c = 0; c < n_channels; ++c) {
		const char *path;

		if (mmap_chan_path != NULL) {
			int ret;

			ret = snprintf(mmap_chan_path, mmap_chan_path_len,
				       "%s.mmap%u", new_path, c);

			if (ret == mmap_chan_path_len) {
				PyErr_SetString(PyExc_SystemError,
						"failed to create mmap path");
				return -1;
			}

			path = mmap_chan_path;
		} else {
			path = open_paths[c];
		}

		if (splat_channel_init_mmap(&frag->channels[c], length, path,
					    frag->temp_mmap))
			return -1;
	}

	if (mmap_chan_path != NULL)
		PyMem_Free(mmap_chan_path);

	return splat_frag_init_common(frag, n_channels, rate, length, name);
}

void splat_frag_free(struct splat_fragment *frag)
{
	unsigned c;
	struct splat_channel *chan;

	for (c = 0, chan = frag->channels; c < frag->n_channels; ++c, ++chan)
		chan->free(chan);

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
	unsigned c;
	struct splat_channel *chan;

	if (length == frag->length)
		return 0;

#if SPLAT_FAST
	length = splat_round4(length);
#endif

	for (c = 0, chan = frag->channels; c < frag->n_channels; ++c, ++chan)
		if (chan->resize(chan, length))
			return -1;

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
		const sf_float_t *src =
			(sf_float_t *)&in->channels[c].data[start];
		sf_float_t *dst =
			(sf_float_t *)&frag->channels[c].data[offset];
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
		const sample_t *src = &incoming->channels[c].data[start];
		sample_t *dst =  &frag->channels[c].data[offset];
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

				frag->channels[c].data[i] +=
					incoming->channels[c].data[in] * a;
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
		sample_t * const chan_data = frag->channels[c].data;
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
		sample_t * const chan_data = frag->channels[c].data;
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
			frag->channels[c].data[i] *= g;
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
				frag->channels[c].data[i] *=
					sig.vectors[c].data[j];
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
			frag->channels[c].data[i] =
				lin2dB(frag->channels[c].data[i]);
	}
}

void splat_frag_dB2lin(struct splat_fragment *frag)
{
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c) {
		size_t i;

		for (i = 0; i < frag->length; ++i)
			frag->channels[c].data[i] =
				dB2lin(frag->channels[c].data[i]);
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
				frag->channels[c].data[i] += offset_float;
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
					frag->channels[c].data[i] += value;
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

static sample_t splat_resample_quad(double x1, const sample_t *y, double x)
{
	const double y0 = y[0];
	const double y1 = y[1];
	const double y2 = y[2];
	double k0, k1, k2;

	k2 = (y2 + y0 - 2.0 * y1) / 2.0;
	k1 = y1 - y0 + k2 * (1.0 - 2.0 * x1);
	k0 = y1 + x1 * (y0 - y1) + k2 * (x1 * x1 - x1);

	return k0 + x * k1 + x * x * k2;
}

int splat_frag_resample_float(struct splat_fragment *frag,
			      const struct splat_fragment *old_frag,
			      unsigned rate, double time_ratio)
{
	const double ratio = time_ratio * (double)rate / frag->rate;
	const size_t max_x0 = old_frag->length - 2;
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
		sample_t *to = frag->channels[c].data;
		size_t i;

		*to++ = old_frag->channels[c].data[0];

		for (i = 1; i < frag->length; ++i) {
			const double x = i / ratio;
			const size_t x0 = minmax(x + 0.5, 1, max_x0);
			const sample_t *y =
				&old_frag->channels[c].data[x0 - 1];

			*to++ = splat_resample_quad(x0, y, x);
		}
	}

	return 0;
}

static int splat_frag_resample_signals(struct splat_fragment *frag,
                                       const struct splat_fragment *old_frag,
                                       unsigned rate, PyObject *ratio)
{
	const double rate_ratio = (double)rate / frag->rate;
	const size_t max_x0 = old_frag->length - 2;
	struct splat_signal sig;
	size_t i;
	double x;
	size_t new_length;
	unsigned c;
	int stop;

	if (splat_signal_init(&sig, frag->length, 0.0, &ratio, 1, rate))
		return -1;

	frag->rate = rate;
	new_length = sig.length;
	stop = 0;

	for (c = 0; c < frag->n_channels; ++c)
		frag->channels[c].data[0] = old_frag->channels[c].data[0];

	i = 1;
	x = 0.0;

	while (!stop && (splat_signal_next(&sig) == SPLAT_SIGNAL_CONTINUE)) {
		size_t j;

		for (j = 0; j < sig.len; ++j, ++i) {
			const double r = sig.vectors[0].data[j] * rate_ratio;
			size_t x0;

			if (r <= 0.0) {
				PyErr_SetString(PyExc_ValueError,
				       "resample time ratio must be positive");
				goto error_free_sig;
			}

			x += 1.0 / r;

			if (i >= frag->length)
				if (splat_frag_resize(frag, frag->length * 1.5))
					goto error_free_sig;

			x0 = minmax(x, 1, max_x0);

			for (c = 0; c < frag->n_channels; ++c) {
				const sample_t *y =
					&old_frag->channels[c].data[x0 - 1];

				frag->channels[c].data[i] =
					splat_resample_quad(x0, y, x);
			}

			if (x0 == max_x0) {
				stop = 1;
				new_length = i;
				break;
			}
		}
	}

	splat_frag_resize(frag, new_length);
	splat_signal_free(&sig);

	return (sig.stat == SPLAT_SIGNAL_ERROR) ? -1 : 0;

error_free_sig:
	splat_signal_free(&sig);

	return -1;
}

static int splat_frag_copy_file_pages(const char *path_from,
				      const char *path_to)
{
	FILE *from;
	FILE *to;
	void *buf;
	int res = -1;

	from = fopen(path_from, "rb");

	if (from == NULL)
		goto out;

	to = fopen(path_to, "wb");

	if (to == NULL)
		goto out_close_from;

	buf = malloc(splat_page_size);

	if (buf == NULL)
		goto out_close_to;

	while (fread(buf, splat_page_size, 1, from) == 1)
		if (fwrite(buf, splat_page_size, 1, to) != 1)
			goto out_free_buf;

	res = 0;

out_free_buf:
	free(buf);
out_close_to:
	fclose(to);
out_close_from:
	fclose(from);
out:
	return res;
}

static int splat_frag_resample_mmap_named(struct splat_fragment *frag,
					  struct splat_fragment *tmp_frag)
{
	const char *paths[SPLAT_MAX_CHANNELS];
	int persist[SPLAT_MAX_CHANNELS];
	unsigned c;
	int res = -1;

	for (c = 0; c < frag->n_channels; ++c)
		paths[c] = NULL;

	for (c = 0; c < frag->n_channels; ++c) {
		persist[c] = frag->channels[c].mmap.persist;
		frag->channels[c].mmap.persist = 1;
		paths[c] = strdup(frag->channels[c].mmap.path);

		if (paths[c] == NULL)
			goto out;
	}

	splat_frag_free(frag);

	for (c = 0; c < tmp_frag->n_channels; ++c) {
		res = unlink(paths[c]);

		if (res)
			goto out;

		res = link(tmp_frag->channels[c].mmap.path, paths[c]);

		if (res && (errno == EXDEV))
			res = splat_frag_copy_file_pages(
				tmp_frag->channels[c].mmap.path, paths[c]);

		if (res)
			goto out;
	}

	res = splat_frag_init_mmap(frag, tmp_frag->n_channels, tmp_frag->rate,
				   tmp_frag->length, tmp_frag->name,
				   NULL, paths);

	if (!res)
		for (c = 0; c < frag->n_channels; ++c)
			frag->channels[c].mmap.persist = persist[c];

out:
	splat_frag_free(tmp_frag);

	for (c = 0; c < frag->n_channels; ++c)
		if (paths[c] != NULL)
			free((char *)paths[c]);

	return res;
}

static int splat_frag_resample_mmap(struct splat_fragment *frag, unsigned rate,
				    PyObject *ratio)
{
	struct splat_fragment tmp_frag;
	int res;

	res = splat_frag_init_mmap(&tmp_frag, frag->n_channels, frag->rate,
				   frag->length, frag->name, NULL, NULL);

	if (res)
		return -1;

	if (PyFloat_Check(ratio))
		res = splat_frag_resample_float(&tmp_frag, frag, rate,
						PyFloat_AsDouble(ratio));
	else
		res = splat_frag_resample_signals(&tmp_frag, frag, rate, ratio);

	if (res)
		return -1;

	if (frag->temp_mmap) {
		splat_frag_free(frag);
		memcpy(frag, &tmp_frag, sizeof(struct splat_fragment));
	} else {
		res = splat_frag_resample_mmap_named(frag, &tmp_frag);
	}

	return res;
}

int splat_frag_resample(struct splat_fragment *frag, unsigned rate,
			PyObject *ratio)
{
	struct splat_fragment old_frag;
	unsigned c;
	int res;

	if (frag->uses_mmap)
		return splat_frag_resample_mmap(frag, rate, ratio);

	if (splat_frag_init(&old_frag, frag->n_channels, frag->rate,
			    frag->length, NULL))
		return -1;

	for (c = 0; c < frag->n_channels; ++c) {
		sample_t *mv;

		mv = old_frag.channels[c].data;
		old_frag.channels[c].data = frag->channels[c].data;
		frag->channels[c].data = mv;
	}

	if (PyFloat_Check(ratio))
		res = splat_frag_resample_float(frag, &old_frag, rate,
						PyFloat_AsDouble(ratio));
	else
		res = splat_frag_resample_signals(frag, &old_frag, rate, ratio);

	splat_frag_free(&old_frag);

	return res;
}

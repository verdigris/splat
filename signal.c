/*
    Splat - signal.c

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

static int splat_signal_func(struct splat_signal *s, struct splat_vector *v)
{
	const double rate = s->rate;
	sample_t *out = v->data;
	size_t i = s->cur;
	size_t j = s->len;

	while (j--) {
		PyObject *ret;

		PyFloat_AS_DOUBLE(s->py_float) = i++ / rate;
		ret = PyObject_Call(v->obj, s->py_args, NULL);

		if (!PyFloat_Check(ret)) {
			PyErr_SetString(PyExc_TypeError,
					"Signal did not return a float");
			Py_DECREF(ret);
			return -1;
		}

		*out++ = PyFloat_AS_DOUBLE(ret);
		Py_DECREF(ret);
	}

	return 0;
}

static int splat_signal_frag(struct splat_signal *s, struct splat_vector *v)
{
	struct splat_fragment *frag = splat_frag_from_obj(v->obj);

	memcpy(v->data, &frag->data[0][s->cur], (s->len * sizeof(sample_t)));

	return 0;
}

static int splat_signal_spline(struct splat_signal *s, struct splat_vector *v)
{
	struct splat_spline *spline = splat_spline_from_obj(v->obj);
	const double rate = s->rate;
	const double k0 = spline->k0;
	sample_t *out = v->data;
	size_t i = s->cur;
	size_t j = s->len;
	PyObject *poly = NULL;
	double end = 0.0;

	while (j--) {
		const double x = i++ / rate;

		if ((x > end) || (poly == NULL)) {
			poly = splat_spline_find_poly(spline->pols, x, &end);

			if (poly == NULL) {
				PyErr_SetString(PyExc_ValueError,
						"Spline polynomial not found");
				return -1;
			}
		}

		*out++ = splat_spline_tuple_value(poly, x) * k0;
	}

	return 0;
}

static int splat_signal_cache(struct splat_signal *s, size_t cur)
{
	size_t i;

	if (cur >= s->length)
		return SPLAT_SIGNAL_STOP;

	if ((cur == s->cur) && s->len)
		return SPLAT_SIGNAL_CONTINUE;

	s->cur = cur;
	s->end = min((s->cur + SPLAT_VECTOR_LEN), s->length);
	s->len = s->end - s->cur;

	for (i = 0; i < s->n_vectors; ++i) {
		struct splat_vector *v = &s->vectors[i];

		if ((v->signal != NULL) && (v->signal(s, v)))
			return SPLAT_SIGNAL_ERROR;
	}

	return SPLAT_SIGNAL_CONTINUE;
}

/* ----------------------------------------------------------------------------
 * Public interface
 */

int splat_signal_init(struct splat_signal *s, size_t length,
		      size_t origin, PyObject **signals,
		      size_t n_signals, unsigned rate)
{
	size_t i;

	s->origin = origin;
	s->length = length + s->origin;
	s->n_vectors = n_signals;
	s->vectors = PyMem_Malloc(n_signals * sizeof(struct splat_vector));
	s->rate = rate;

	if (s->vectors == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	s->py_float = PyFloat_FromDouble(0);

	if (s->py_float == NULL) {
		PyErr_SetString(PyExc_AssertionError,
				"Failed to create float object");
		return -1;
	}

	s->py_args = PyTuple_New(1);

	if (s->py_args == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	PyTuple_SET_ITEM(s->py_args, 0, s->py_float);

	for (i = 0; i < n_signals; ++i) {
		struct splat_vector *v = &s->vectors[i];
		PyObject *signal = signals[i];
		struct splat_spline *spline = splat_spline_from_obj(signal);
		struct splat_fragment *frag = splat_frag_from_obj(signal);

		if (PyFloat_Check(signal)) {
			const sample_t value = PyFloat_AS_DOUBLE(signal);
			size_t j;

			for (j = 0; j < SPLAT_VECTOR_LEN; ++j)
				v->data[j] = value;

			v->signal = NULL;
		} else if (PyCallable_Check(signal)) {
			v->signal = splat_signal_func;
		} else if (frag != NULL) {
			if (frag->n_channels != 1) {
				PyErr_SetString(PyExc_ValueError,
				"Fragment signal must have only 1 channel");
				return -1;
			}

			if (s->length > frag->length) {
				PyErr_SetString(PyExc_ValueError,
				"Fragment signal length too short");
				return -1;
			}

			v->signal = splat_signal_frag;
		} else if (spline != NULL) {
			size_t spline_length = spline->end * rate;

			if (s->length > spline_length) {
				PyErr_SetString(PyExc_ValueError,
				"Spline signal length too short");
				return -1;
			}

			v->signal = splat_signal_spline;
		} else {
			PyErr_SetString(PyExc_TypeError,
					"unsupported signal type");
			return -1;
		}

		v->obj = signal;
	}

	s->cur = s->origin;
	s->end = 0;
	s->len = 0;
	s->stat = SPLAT_SIGNAL_CONTINUE;

	return 0;
}

void splat_signal_free(struct splat_signal *s)
{
	Py_DECREF(s->py_float);
	Py_DECREF(s->py_args);
	PyMem_Free(s->vectors);
}

int splat_signal_next(struct splat_signal *s)
{
	s->stat = splat_signal_cache(s, (s->cur + s->len));

	return s->stat;
}

ssize_t splat_signal_get(struct splat_signal *s, size_t n)
{
	div_t co; /* cursor, offset */
	size_t cur;

	if (n >= s->length)
		return -1;

	co = div(n, SPLAT_VECTOR_LEN);
	cur = co.quot * SPLAT_VECTOR_LEN;
	s->stat = splat_signal_cache(s, cur);

	if (s->stat != SPLAT_SIGNAL_CONTINUE)
		return -1;

	return co.rem;
}

PyObject *splat_signal_tuple(struct splat_signal *s, size_t offset)
{
	PyObject *sig_tuple;
	size_t i;

	sig_tuple = PyTuple_New(s->n_vectors);

	if (sig_tuple == NULL)
		return NULL;

	for (i = 0; i < s->n_vectors; ++i) {
		PyObject *val = PyFloat_FromDouble(s->vectors[i].data[offset]);

		if (val == NULL) {
			Py_DECREF(sig_tuple);
			return NULL;
		}

		PyTuple_SET_ITEM(sig_tuple, i, val);
	}

	return sig_tuple;
}

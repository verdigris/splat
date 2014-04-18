/*
    Splat - _splat.c

    Copyright (C) 2012, 2013, 2014 Guillaume Tucker <guillaume@mangoz.org>

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

#include <Python.h>
#include <math.h>

/* Set to 1 to use 4xfloat vectors (SSE) */
#define USE_V4SF 0

#define BASE_TYPE_FLAGS (Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE)
#define MAX_CHANNELS 16

#define lin2dB(level) (20 * log10(level))
#define dB2lin(dB) (pow10((dB) / 20))

/* Use PP_NARG to get the number of arguments in __VA_ARGS__ */
#define PP_NARG(...) \
         PP_NARG_(__VA_ARGS__,PP_RSEQ_N())
#define PP_NARG_(...) \
         PP_ARG_N(__VA_ARGS__)
#define PP_ARG_N( \
          _1, _2, _3, _4, _5, _6, _7, _8, _9,_10, \
         _11,_12,_13,_14,_15,_16,_17,_18,_19,_20, \
         _21,_22,_23,_24,_25,_26,_27,_28,_29,_30, \
         _31,_32,_33,_34,_35,_36,_37,_38,_39,_40, \
         _41,_42,_43,_44,_45,_46,_47,_48,_49,_50, \
         _51,_52,_53,_54,_55,_56,_57,_58,_59,_60, \
         _61,_62,_63,N,...) N
#define PP_RSEQ_N() \
         63,62,61,60,                   \
         59,58,57,56,55,54,53,52,51,50, \
         49,48,47,46,45,44,43,42,41,40, \
         39,38,37,36,35,34,33,32,31,30, \
         29,28,27,26,25,24,23,22,21,20, \
         19,18,17,16,15,14,13,12,11,10, \
         9,8,7,6,5,4,3,2,1,0

#ifndef min
# define min(_a, _b) (((_a) < (_b)) ? (_a) : (_b))
#endif

#ifndef max
# define max(_a, _b) (((_a) > (_b)) ? (_a) : (_b))
#endif

#ifndef minmax
# define minmax(_val, _min, _max) \
	((_val) < (_min) ? (_min) : ((_val) > (_max) ? (_max) : (_val)))
#endif

/* Type of each sample */
#if USE_V4SF /* for speed and smaller memory footprint */
typedef float sample_t;
typedef float v4sf __attribute__ ((vector_size(16)));
#else /* for precision */
typedef double sample_t;
#endif

/* ----------------------------------------------------------------------------
 * Module constants
 */

/* Initial ratio for sound sources (square, triangle) */
static PyObject *splat_init_source_ratio;

/* ----------------------------------------------------------------------------
 * Fragment class interface
 */

struct Fragment_object {
	PyObject_HEAD;
	int init;
	unsigned n_channels;
	unsigned rate;
	Py_ssize_t length;
	sample_t *data[MAX_CHANNELS];
};
typedef struct Fragment_object Fragment;

static PyTypeObject splat_FragmentType;


/* ----------------------------------------------------------------------------
 * Signal vector interface
 */

#define SIGNAL_VECTOR_BITS 8
#define SIGNAL_VECTOR_LEN (1 << SIGNAL_VECTOR_BITS)
#define SIGNAL_VECTOR_MASK (SIGNAL_VECTOR_LEN - 1)

struct splat_signal;

struct signal_vector {
	sample_t data[SIGNAL_VECTOR_LEN];
	PyObject *obj;
	int (*signal)(struct splat_signal *s, struct signal_vector *v);
};

enum signal_ret {
	SIGNAL_VECTOR_CONTINUE = 0,
	SIGNAL_VECTOR_STOP,
};

struct splat_signal {
	Fragment *frag;
	size_t n_vectors;
	struct signal_vector *vectors;
	PyObject *py_float;
	PyObject *py_args;
	size_t cur;
	size_t end;
	size_t len;
	int error;
};

static int splat_signal_init(struct splat_signal *s, Fragment *frag,
			     PyObject **signals, size_t n_signals);
static void splat_signal_free(struct splat_signal *s);
static int splat_signal_next(struct splat_signal *s);
static ssize_t splat_signal_get(struct splat_signal *s, size_t n);
static PyObject *splat_signal_tuple(struct splat_signal *s, size_t offset);

/* -- Signal class -- */

struct Signal_object {
	PyObject_HEAD;
	int init;
	PyObject **signals;
	struct splat_signal sig;
	enum signal_ret stat;
	unsigned offset;
};
typedef struct Signal_object Signal;

static PyTypeObject splat_SignalType;

static void Signal_free_signals(Signal *self, size_t n)
{
	size_t i;

	for (i = 0; i < n; ++i)
		Py_DECREF(self->signals[i]);

	PyMem_Free(self->signals);
}

static void Signal_dealloc(Signal *self)
{
	if (self->init) {
		splat_signal_free(&self->sig);
		Signal_free_signals(self, self->sig.n_vectors);
		Py_DECREF(self->sig.frag);
		self->init = 0;
	}

	self->ob_type->tp_free((PyObject *)self);
}

static int Signal_init(Signal *self, PyObject *args)
{
	Fragment *frag;
	PyObject *sig_obj;

	size_t n_signals;
	size_t i;

	if (!PyArg_ParseTuple(args, "O!O", &splat_FragmentType, &frag,
			      &sig_obj))
		return -1;

	if (PyTuple_Check(sig_obj))
		n_signals = PyTuple_GET_SIZE(sig_obj);
	else
		n_signals = 1;

	self->signals = PyMem_Malloc(sizeof(PyObject *) * n_signals);

	if (self->signals == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	if (PyTuple_Check(sig_obj)) {
		for (i = 0; i < n_signals; ++i)
			self->signals[i] = PyTuple_GET_ITEM(sig_obj, i);
	} else {
		self->signals[0] = sig_obj;
	}

	for (i = 0; i < n_signals; ++i)
		Py_INCREF(self->signals[i]);

	if (splat_signal_init(&self->sig, frag, self->signals, n_signals)) {
		Signal_free_signals(self, n_signals);
		return -1;
	}

	Py_INCREF(self->sig.frag);
	self->offset = 0;
	self->stat = SIGNAL_VECTOR_CONTINUE;
	self->init = 1;

	return 0;
}

static PyObject *Signal_new(PyTypeObject *type, PyObject *args, PyObject *kw)
{
	Signal *self;

	self = (Signal *)type->tp_alloc(type, 0);

	if (self == NULL)
		return PyErr_NoMemory();

	self->init = 0;

	return (PyObject *)self;
}

/* Signal sequence interface */

static Py_ssize_t SignalObj_sq_length(Signal *self)
{
	return self->sig.frag->length;
}

static PyObject *SignalObj_sq_item(Signal *self, Py_ssize_t i)
{
	ssize_t offset;

	offset = splat_signal_get(&self->sig, i);

	if (offset < 0)
		return NULL;

	return splat_signal_tuple(&self->sig, offset);
}

static PySequenceMethods Signal_as_sequence = {
	(lenfunc)SignalObj_sq_length, /* sq_length (lenfunc) */
	NULL, /* sq_concat (binaryfunc) */
	NULL, /* sq_repeat (ssizeargfunc) */
	(ssizeargfunc)SignalObj_sq_item, /* sq_item (ssizeargfunc) */
	NULL, /* sq_slice (ssizessizeargfunc) */
	NULL, /* sq_ass_item (ssizeobjargproc) */
	NULL, /* sq_ass_slice (ssizessizeobjargproc) */
	NULL, /* sq_contains (objobjproc) */
	NULL, /* sq_inplace_concat (binaryfunc) */
	NULL, /* sq_inplace_repeat (ssizeargfunc) */
};

static PyTypeObject splat_SignalType = {
	PyObject_HEAD_INIT(NULL)
	0,                                 /* ob_size */
	"_splat.Signal",                   /* tp_name */
	sizeof(Signal),                    /* tp_basicsize */
	0,                                 /* tp_itemsize */
	(destructor)Signal_dealloc,        /* tp_dealloc */
	0,                                 /* tp_print */
	0,                                 /* tp_getattr */
	0,                                 /* tp_setattr */
	0,                                 /* tp_compare */
	0,                                 /* tp_repr */
	0,                                 /* tp_as_number */
	&Signal_as_sequence,               /* tp_as_sequence */
	0,                                 /* tp_as_mapping */
	0,                                 /* tp_hash  */
	0,                                 /* tp_call */
	0,                                 /* tp_str */
	0,                                 /* tp_getattro */
	0,                                 /* tp_setattro */
	0,                                 /* tp_as_buffer */
	BASE_TYPE_FLAGS,                   /* tp_flags */
	"Signal object",                   /* tp_doc */
	0,                                 /* tp_traverse */
	0,                                 /* tp_clear */
	0,                                 /* tp_richcompare */
	0,                                 /* tp_weaklistoffset */
	0,                                 /* tp_iter */
	0,                                 /* tp_iternext */
	0,                                 /* tp_methods */
	0,                                 /* tp_members */
	0,                                 /* tp_getset */
	0,                                 /* tp_base */
	0,                                 /* tp_dict */
	0,                                 /* tp_descr_get */
	0,                                 /* tp_descr_set */
	0,                                 /* tp_dictoffset */
	(initproc)Signal_init,             /* tp_init */
	0,                                 /* tp_alloc */
	Signal_new,                        /* tp_new */
};

/* -- Signal internal functions -- */

static int splat_signal_func(struct splat_signal *s, struct signal_vector *v)
{
	const double rate = s->frag->rate;
	size_t i;

	for (i = s->cur; i < s->end; ++i) {
		PyObject *ret;

		PyFloat_AS_DOUBLE(s->py_float) = i / rate;
		ret = PyObject_Call(v->obj, s->py_args, NULL);

		if (!PyFloat_Check(ret))
			return -1;

		v->data[i & SIGNAL_VECTOR_MASK] = PyFloat_AS_DOUBLE(ret);
	}

	return 0;
}

static int splat_signal_frag(struct splat_signal *s, struct signal_vector *v)
{
	Fragment *frag = (Fragment *)v->obj;
	size_t i;

	for (i = s->cur; i < s->end; ++i)
		v->data[i & SIGNAL_VECTOR_MASK] = frag->data[0][i];

	return 0;
}

static int splat_signal_init(struct splat_signal *s, Fragment *frag,
			     PyObject **signals, size_t n_signals)
{
	size_t i;

	s->frag = frag;
	s->n_vectors = n_signals;
	s->vectors = PyMem_Malloc(n_signals * sizeof(struct signal_vector));

	if (s->vectors == NULL) {
		PyErr_NoMemory();
		return -1;
	}

	s->py_float = PyFloat_FromDouble(0);

	if (s->py_float == NULL)
		return -1;

	s->py_args = PyTuple_New(1);

	if (s->py_args == NULL)
		return -1;

	PyTuple_SET_ITEM(s->py_args, 0, s->py_float);

	for (i = 0; i < n_signals; ++i) {
		struct signal_vector *v = &s->vectors[i];
		PyObject *signal = signals[i];

		if (PyFloat_Check(signal)) {
			const sample_t value = PyFloat_AS_DOUBLE(signal);
			size_t j;

			for (j = 0; j < SIGNAL_VECTOR_LEN; ++j)
				v->data[j] = value;

			v->signal = NULL;
		} else if (PyCallable_Check(signal)) {
			v->signal = splat_signal_func;
		} else if (PyObject_TypeCheck(signal, &splat_FragmentType)) {
			Fragment *sig_frag = (Fragment *)signal;

			if (sig_frag->length != frag->length) {
				PyErr_SetString(PyExc_ValueError,
						"fragment size mismatch");
				return -1;
			}

			v->signal = splat_signal_frag;
		} else {
			PyErr_SetString(PyExc_TypeError,
					"unsupported signal type");
			return -1;
		}

		v->obj = signal;
	}

	s->cur = 0;
	s->end = 0;
	s->len = 0;
	s->error = 0;

	return 0;
}

static void splat_signal_free(struct splat_signal *s)
{
	Py_DECREF(s->py_args);
	PyMem_Free(s->vectors);
}

static int splat_signal_cache(struct splat_signal *s, size_t cur)
{
	size_t i;

	if (cur >= s->frag->length)
		return SIGNAL_VECTOR_STOP;

	if ((cur == s->cur) && s->len)
		return SIGNAL_VECTOR_CONTINUE;

	s->cur = cur;
	s->end = min((s->cur + SIGNAL_VECTOR_LEN), s->frag->length);
	s->len = s->end - s->cur;

	for (i = 0; i < s->n_vectors; ++i) {
		struct signal_vector *v = &s->vectors[i];

		if ((v->signal != NULL) && (v->signal(s, v)))
			return SIGNAL_VECTOR_STOP;
	}

	return SIGNAL_VECTOR_CONTINUE;
}

static int splat_signal_next(struct splat_signal *s)
{
	return splat_signal_cache(s, (s->cur + s->len));
}

static ssize_t splat_signal_get(struct splat_signal *s, size_t n)
{
	div_t co; /* cursor, offset */
	size_t cur;

	if (n >= s->frag->length)
		return -1;

	co = div(n, SIGNAL_VECTOR_LEN);
	cur = co.quot * SIGNAL_VECTOR_LEN;

	if (splat_signal_cache(s, cur) == SIGNAL_VECTOR_STOP)
		return -1;

	return co.rem;
}

static PyObject *splat_signal_tuple(struct splat_signal *s, size_t offset)
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


/* ----------------------------------------------------------------------------
 * Fragment class
 */

/* -- internal functions -- */

struct splat_levels_helper {
	unsigned n;
	PyObject *obj[MAX_CHANNELS];
	double fl[MAX_CHANNELS]; /* levels converted to linear scale */
	int all_floats;
};

static int frag_get_levels(Fragment *frag, struct splat_levels_helper *levels,
			   PyObject *levels_obj);
static int frag_resize(Fragment *self, size_t length);
#define frag_grow(_frag, _length)					\
	(((_length) <= (_frag)->length) ? 0 : frag_resize((_frag), (_length)))

/* -- Fragment functions -- */

static void Fragment_dealloc(Fragment *self)
{
	if (self->init) {
		unsigned i;

		for (i = 0; i < self->n_channels; ++i)
			PyMem_Free(self->data[i]);

		self->init = 0;
	}

	self->ob_type->tp_free((PyObject *)self);
}

static int Fragment_init(Fragment *self, PyObject *args, PyObject *kw)
{
	static char *kwlist[] = { "channels", "rate", "duration", NULL };
	unsigned n_channels = 2;
	unsigned rate = 48000;
	double duration = 0.0;

	unsigned i;
	size_t length;
	size_t data_size;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "|IId", kwlist,
					 &n_channels, &rate, &duration))
		return -1;

	if (duration < 0.0) {
		PyErr_SetString(PyExc_ValueError, "negative duration");
		return -1;
	}

	if (n_channels > MAX_CHANNELS) {
		PyErr_SetString(PyExc_ValueError,
				"exceeding maximum number of channels");
		return -1;
	}

	length = duration * rate;
	data_size = length * sizeof(sample_t);

	for (i = 0; i < n_channels; ++i) {
		if (!data_size) {
			self->data[i] = NULL;
		} else {
			self->data[i] = PyMem_Malloc(data_size);

			if (self->data[i] == NULL) {
				PyErr_NoMemory();
				return -1;
			}

			memset(self->data[i], 0, data_size);
		}
	}

	self->n_channels = n_channels;
	self->rate = rate;
	self->length = length;
	self->init = 1;

	return 0;
}

static PyObject *Fragment_new(PyTypeObject *type, PyObject *args, PyObject *kw)
{
	Fragment *self;

	self = (Fragment *)type->tp_alloc(type, 0);

	if (self == NULL)
		return PyErr_NoMemory();

	self->init = 0;

	return (PyObject *)self;
}

/* Fragment sequence interface */

static Py_ssize_t Fragment_sq_length(Fragment *self)
{
	return self->length;
}

static PyObject *Fragment_sq_item(Fragment *self, Py_ssize_t i)
{
	PyObject *sample;
	unsigned c;

	if ((i < 0) || (i >= self->length)) {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	sample = PyTuple_New(self->n_channels);

	if (sample == NULL)
		return NULL;

	for (c = 0; c < self->n_channels; ++c) {
		const sample_t s = self->data[c][i];
		PyTuple_SetItem(sample, c, PyFloat_FromDouble(s));
	}

	return sample;
}

static int Fragment_sq_ass_item(Fragment *self, Py_ssize_t i, PyObject *v)
{
	unsigned c;

	if (!PyTuple_CheckExact(v)) {
		PyErr_SetString(PyExc_TypeError, "item must be a tuple");
		return -1;
	}

	if (PyTuple_GET_SIZE(v) != self->n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return -1;
	}

	if ((i < 0) || (i >= self->length)) {
		PyErr_SetString(PyExc_IndexError, "set index error");
		return -1;
	}

	for (c = 0; c < self->n_channels; ++c) {
		PyObject *s = PyTuple_GET_ITEM(v, c);

		if (!PyFloat_CheckExact(s)) {
			PyErr_SetString(PyExc_TypeError,
					"item must contain floats");
			return -1;
		}

		self->data[c][i] = (sample_t)PyFloat_AS_DOUBLE(s);
	}

	return 0;
}

static PySequenceMethods Fragment_as_sequence = {
	(lenfunc)Fragment_sq_length, /* sq_length */
	(binaryfunc)0, /* sq_concat */
	(ssizeargfunc)0, /* sq_repeat */
	(ssizeargfunc)Fragment_sq_item, /* sq_item */
	(ssizessizeargfunc)0, /* sq_slice */
	(ssizeobjargproc)Fragment_sq_ass_item, /* sq_ass_item */
	(ssizessizeobjargproc)0, /* sq_ass_slice */
	(objobjproc)0, /* sq_contains */
	(binaryfunc)0, /* sq_inplace_concat */
	(ssizeargfunc)0, /* sq_inplace_repeat */
};

/* Fragment getsetters */

PyDoc_STRVAR(sample_rate_doc, "Get the sample rate in Hz.");

static PyObject *Fragment_get_sample_rate(Fragment *self, void *_)
{
	return Py_BuildValue("I", self->rate);
}

PyDoc_STRVAR(duration_doc, "Get the fragment duration in seconds.");

static PyObject *Fragment_get_duration(Fragment *self, void *_)
{
	return Py_BuildValue("f", (double)self->length / self->rate);
}

PyDoc_STRVAR(channels_doc, "Get the number of channels.");

static PyObject *Fragment_get_channels(Fragment *self, void *_)
{
	return Py_BuildValue("I", self->n_channels);
}

static PyGetSetDef Fragment_getsetters[] = {
	{ "sample_rate", (getter)Fragment_get_sample_rate, NULL,
	  sample_rate_doc },
	{ "duration", (getter)Fragment_get_duration, NULL, duration_doc },
	{ "channels", (getter)Fragment_get_channels, NULL, channels_doc },
	{ NULL }
};

/* Fragment methods */

PyDoc_STRVAR(Fragment_mix_doc,
"mix(fragment, offset=0.0, start=0.0, levels=None, length=None)\n"
"\n"
"Mix the given other ``fragment`` data into this instance.\n"
"\n"
"This is achieved by simply adding up the samples of an incoming fragment "
"to this fragment.  The ``offset``, ``start`` and ``length`` values in "
"seconds can be used to alter the mixing times.  The incoming fragment can "
"start being mixed with an ``offset`` into this fragment, the beginning can "
"be skipped until the given ``start`` time, and the ``length`` to be mixed "
"can be manually limited.  These values will be automatically adjusted to "
"remain within the available data incoming.  The length of this fragment will "
"be automatically increased if necessary to receive the incoming data.\n");

static PyObject *Fragment_mix(Fragment *self, PyObject *args, PyObject *kw)
{
	static char *kwlist[] = {
		"frag", "offset", "start", "levels", "length", NULL };
	Fragment *frag;
	double offset = 0.0;
	double start = 0.0;
	PyObject *length_obj = NULL;
	PyObject *levels_obj = NULL;

	struct splat_levels_helper levels;
	ssize_t length;
	ssize_t offset_sample;
	ssize_t start_sample;
	size_t total_length;
	unsigned c;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|ddOd", kwlist,
					 &splat_FragmentType, &frag, &offset,
					 &start, &levels_obj, &length_obj))
		return NULL;

	if (frag->n_channels != self->n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return NULL;
	}

	if (frag->rate != self->rate) {
		PyErr_SetString(PyExc_ValueError, "sample rate mismatch");
		return NULL;
	}

	if (levels_obj != NULL)
		if (frag_get_levels(self, &levels, levels_obj))
			return NULL;

	if (length_obj != NULL) {
		if (!PyFloat_Check(length_obj)) {
			PyErr_SetString(PyExc_ValueError,
					"length must be float");
			return NULL;
		}

		length = PyFloat_AS_DOUBLE(length_obj) * self->rate;
	} else {
		length = frag->length;
	}

	offset_sample = offset * self->rate;
	offset_sample = max(offset_sample, 0);
	start_sample = start * self->rate;
	start_sample = minmax(start_sample, 0, frag->length);
	length = minmax(length, 0, (frag->length - start_sample));
	total_length = offset_sample + length;

	if (frag_grow(self, total_length))
		return NULL;

	for (c = 0; c < self->n_channels; ++c) {
		const sample_t *src = &frag->data[c][start_sample];
		sample_t *dst =  &self->data[c][offset_sample];
		Py_ssize_t i = length;

		if (levels_obj == NULL) {

			while (i--)
				*dst++ += *src++;
		} else {
			const double g = levels.fl[c];

			while (i--)
				*dst++ += g * (*src++);
		}
	}

	Py_RETURN_NONE;
}

static PyObject *Fragment_import_bytes(Fragment *self, PyObject *args)
{
	PyObject *bytes_obj;
	int start;
	unsigned sample_width;
	unsigned sample_rate;
	unsigned n_channels;

	const char *bytes;
	Py_ssize_t n_bytes;
	unsigned sample_bits;
	unsigned bytes_per_sample;
	unsigned n_samples;
	unsigned end;
	unsigned c;
	int32_t neg;
	int32_t neg_mask;
	double scale;

	if (!PyArg_ParseTuple(args, "O!iIII", &PyByteArray_Type, &bytes_obj,
			      &start, &sample_width, &sample_rate,
			      &n_channels))
		return NULL;

	if ((sample_width != 2) && (sample_width != 3)) {
		PyErr_SetString(PyExc_ValueError, "unsupported sample width");
		return NULL;
	}

	if (n_channels != self->n_channels) {
		PyErr_SetString(PyExc_ValueError, "wrong number of channels");
		return NULL;
	}

	if (sample_rate != self->rate) {
		PyErr_SetString(PyExc_ValueError, "wrong sample rate");
		return NULL;
	}

	n_bytes = PyByteArray_Size(bytes_obj);
	bytes_per_sample = n_channels * sample_width;

	if (n_bytes % bytes_per_sample) {
		PyErr_SetString(PyExc_ValueError, "invalid buffer length");
		return NULL;
	}

	bytes = PyByteArray_AsString(bytes_obj);
	n_samples = n_bytes / bytes_per_sample;
	end = start + n_samples;

	if (frag_grow(self, end))
		return NULL;

	sample_bits = 8 * sample_width;
	scale = pow(2, sample_bits) / 2;
	neg = 1 << (sample_bits - 1);
	neg_mask = 0xFFFFFFFF - (1 << sample_bits) + 1;

	for (c = 0; c < self->n_channels; ++c) {
		const char *in = bytes + (sample_width * c);
		sample_t *out = &self->data[c][start];
		unsigned s;

		if (sample_width == 2) {
			for (s = start; s < end; ++s) {
				*out++ = *(int16_t *)in / scale;
				in += bytes_per_sample;
			}
		} else {
			for (s = start; s < end; ++s) {
				const uint8_t *b = (const uint8_t *)in;
				int32_t sample = 0;
				int n;

				for (n = 0; n < sample_width; ++n)
					sample += (*b++) << (n * 8);

				if (sample & neg)
					sample |= neg_mask;

				*out++ = sample / scale;
				in += bytes_per_sample;
			}
		}
	}

	Py_RETURN_NONE;
}

static PyObject *Fragment_as_bytes(Fragment *self, PyObject *args)
{
	unsigned sample_width;

	PyObject *bytes_obj;
	Py_ssize_t bytes_size;
	const sample_t *it[MAX_CHANNELS];
	unsigned i;
	unsigned c;
	char *out;

	if (!PyArg_ParseTuple(args, "I", &sample_width))
		return NULL;

	if ((sample_width != 1) && (sample_width != 2)) {
		PyErr_SetString(PyExc_ValueError, "unsupported sample width");
		return NULL;
	}

	bytes_obj = PyByteArray_FromStringAndSize("", 0);

	if (bytes_obj == NULL)
		return PyErr_NoMemory();

	bytes_size = self->length * self->n_channels * sample_width;

	if (PyByteArray_Resize(bytes_obj, bytes_size))
		return PyErr_NoMemory();

	for (c = 0; c < self->n_channels; ++c)
		it[c] = self->data[c];

	out = PyByteArray_AS_STRING(bytes_obj);

	if (sample_width == 2) {
		for (i = 0; i < self->length; ++i) {
			for (c = 0; c < self->n_channels; ++c) {
				const sample_t z = *(it[c]++);
				int16_t s;

				if (z < -1.0)
					s = -32767;
				else if (z > 1.0)
					s = 32767;
				else
					s = z * 32767;

				*out++ = s & 0xFF;
				*out++ = (s >> 8) & 0xFF;
			}
		}
	} else if (sample_width == 1) {
		for (i = 0; i < self->length; ++i) {
			for (c = 0; c < self->n_channels; ++c) {
				const sample_t z = *(it[c]++);
				int8_t s;

				if (z < -1.0)
					s = 0;
				else if (z > 1.0)
					s = 255;
				else
					s = (z * 127) + 128;

				*out++ = s & 0xFF;
			}
		}
	}

	return bytes_obj;
}

PyDoc_STRVAR(Fragment_normalize_doc,
"normalize(level=-0.05, zero=True)\n"
"\n"
"Normalize the amplitude.\n"
"\n"
"The ``level`` value in dB is the resulting maximum amplitude after "
"normalization.  The default value of -0.05 dB is the maximum level while "
"ensuring no clipping occurs, even when saving with 8-bit sample resolution. "
"The same gain is applied to all channels, so the relative difference in "
"levels between channels is preserved.\n"
"\n"
"When ``zero`` is ``True``, the average value is substracted from all the "
"fragment prior to amplification to avoid any offset and achieve maximum "
"amplitude.  With some imbalanced transitory signals, it may be better to not "
"remove the average value as this may have the undesirable effect of adding "
"some offset instead.\n");

static PyObject *Fragment_normalize(Fragment *self, PyObject *args)
{
	double level = -0.05;
	PyObject *zero = NULL;

	int do_zero;
	double average[MAX_CHANNELS];
	unsigned c;
	double peak;
	double gain;

	if (!PyArg_ParseTuple(args, "|dO!", &level, &PyBool_Type, &zero))
		return NULL;

	if (self->n_channels > MAX_CHANNELS) {
		PyErr_SetString(PyExc_ValueError, "too many channels");
		return NULL;
	}

	level = dB2lin(level);
	peak = 0.0;
	do_zero = ((zero == NULL) || (zero == Py_True)) ? 1 : 0;

	for (c = 0; c < self->n_channels; ++c) {
		sample_t * const chan_data = self->data[c];
		const sample_t * const end = &chan_data[self->length];
		const sample_t *it;
		double avg = 0.0;
		double neg = 1.0;
		double pos = -1.0;
		double chan_peak;

		for (it = self->data[c]; it != end; ++it) {
			if (do_zero)
				avg += *it / self->length;

			if (*it > pos)
				pos = *it;

			if (*it < neg)
				neg = *it;
		}

		average[c] = avg;

		if (do_zero) {
			neg -= avg;
			pos -= avg;
		}

		neg = fabsf(neg);
		pos = fabsf(pos);
		chan_peak = (neg > pos) ? neg : pos;

		if (chan_peak > peak)
			peak = chan_peak;
	}

	gain = level / peak;

	for (c = 0; c < self->n_channels; ++c) {
		const double chan_avg = average[c];
		sample_t * const chan_data = self->data[c];
		const sample_t * const end = &chan_data[self->length];
		sample_t *it;

		for (it = chan_data; it != end; ++it) {
			*it -= chan_avg;
			*it *= gain;
		}
	}

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_amp_doc,
"amp(gain)\n"
"\n"
"Amplify the fragment by the given ``gain`` in dB which can either be a "
"single floating point value to apply to all channels or a tuple with a value "
"for each individual channel.\n");

static PyObject *Fragment_amp(Fragment *self, PyObject *args)
{
	PyObject *gain_obj;

	struct splat_levels_helper gains;
	unsigned c;

	if (!PyArg_ParseTuple(args, "O", &gain_obj))
		return NULL;

	if (frag_get_levels(self, &gains, gain_obj))
		return NULL;

	for (c = 0; c < self->n_channels; ++c) {
		const double g = gains.fl[c];
		size_t i;

		for (i = 0; i < self->length; ++i)
			self->data[c][i] *= g;
	}

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_offset_doc,
"offset(value)\n"
"\n"
"Add an offset to the data already in the fragment.  This is especially "
"useful when generating a modulation fragment.\n");

static PyObject *Fragment_offset(Fragment *self, PyObject *args)
{
	PyObject *offset;

	unsigned c;
	size_t i;

	if (!PyArg_ParseTuple(args, "O", &offset))
		return NULL;

	if (PyFloat_Check(offset)) {
		const double offset_float = PyFloat_AS_DOUBLE(offset);

		for (c = 0; c < self->n_channels; ++c) {
			for (i = 0; i < self->length; ++i)
				self->data[c][i] += offset_float;
		}
	} else {
		struct splat_signal sig;

		if (splat_signal_init(&sig, self, &offset, 1))
			return NULL;

		while (splat_signal_next(&sig) == SIGNAL_VECTOR_CONTINUE) {
			size_t j;

			for (i = sig.cur, j = 0; i < sig.end; ++i, ++j) {
				const double value = sig.vectors[0].data[j];

				for (c = 0; c < self->n_channels; ++c)
					self->data[c][i] += value;
			}
		}
	}

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_resize_doc,
"resize(duration)\n"
"\n"
"Resize the fragment to the given ``duration`` in seconds.  If the fragment "
"grows, silence is added to the end.  When shrinking, the end of the fragment "
"is lost.\n");

static PyObject *Fragment_resize(Fragment *self, PyObject *args)
{
	double duration;

	ssize_t length;

	if (!PyArg_ParseTuple(args, "d", &duration))
		return NULL;

	length = duration * self->rate;

	if (frag_resize(self, length))
		return NULL;

	Py_RETURN_NONE;
}

static PyMethodDef Fragment_methods[] = {
	{ "mix", (PyCFunction)Fragment_mix, METH_KEYWORDS,
	  Fragment_mix_doc },
	{ "import_bytes", (PyCFunction)Fragment_import_bytes, METH_VARARGS,
	  "Import data as raw bytes" },
	{ "as_bytes", (PyCFunction)Fragment_as_bytes, METH_VARARGS,
	  "Make a byte buffer with the data" },
	{ "normalize", (PyCFunction)Fragment_normalize, METH_VARARGS,
	  Fragment_normalize_doc },
	{ "amp", (PyCFunction)Fragment_amp, METH_VARARGS,
	  Fragment_amp_doc },
	{ "offset", (PyCFunction)Fragment_offset, METH_VARARGS,
	  Fragment_offset_doc },
	{ "resize", (PyCFunction)Fragment_resize, METH_VARARGS,
	  Fragment_resize_doc },
	{ NULL }
};

static PyTypeObject splat_FragmentType = {
	PyObject_HEAD_INIT(NULL)
	0,                                 /* ob_size */
	"_splat.Fragment",                 /* tp_name */
	sizeof(Fragment),                  /* tp_basicsize */
	0,                                 /* tp_itemsize */
	(destructor)Fragment_dealloc,      /* tp_dealloc */
	0,                                 /* tp_print */
	0,                                 /* tp_getattr */
	0,                                 /* tp_setattr */
	0,                                 /* tp_compare */
	0,                                 /* tp_repr */
	0,                                 /* tp_as_number */
	&Fragment_as_sequence,             /* tp_as_sequence */
	0,                                 /* tp_as_mapping */
	0,                                 /* tp_hash  */
	0,                                 /* tp_call */
	0,                                 /* tp_str */
	0,                                 /* tp_getattro */
	0,                                 /* tp_setattro */
	0,                                 /* tp_as_buffer */
	BASE_TYPE_FLAGS,                   /* tp_flags */
	"Fragment of audio data",          /* tp_doc */
	0,                                 /* tp_traverse */
	0,                                 /* tp_clear */
	0,                                 /* tp_richcompare */
	0,                                 /* tp_weaklistoffset */
	0,                                 /* tp_iter */
	0,                                 /* tp_iternext */
	Fragment_methods,                  /* tp_methods */
	0,                                 /* tp_members */
	Fragment_getsetters,               /* tp_getset */
	0,                                 /* tp_base */
	0,                                 /* tp_dict */
	0,                                 /* tp_descr_get */
	0,                                 /* tp_descr_set */
	0,                                 /* tp_dictoffset */
	(initproc)Fragment_init,           /* tp_init */
	0,                                 /* tp_alloc */
	Fragment_new,                      /* tp_new */
};

/* -- Fragment internal functions -- */

static int frag_get_levels(Fragment *frag, struct splat_levels_helper *levels,
			   PyObject *levels_obj)
{
	unsigned c;

	if (PyFloat_Check(levels_obj)) {
		const double gain_lin = dB2lin(PyFloat_AsDouble(levels_obj));

		for (c = 0; c < frag->n_channels; ++c) {
			levels->obj[c] = levels_obj;
			levels->fl[c] = gain_lin;
			levels->all_floats = 1;
			levels->n = frag->n_channels;
		}
	} else if (PyTuple_Check(levels_obj)) {
		const Py_ssize_t n_channels = PyTuple_GET_SIZE(levels_obj);

		if (n_channels > MAX_CHANNELS) {
			PyErr_SetString(PyExc_ValueError, "too many channels");
			return -1;
		}

		if (n_channels != frag->n_channels) {
			PyErr_SetString(PyExc_ValueError,
					"channels number mismatch");
			return -1;
		}

		levels->n = n_channels;
		levels->all_floats = 1;

		for (c = 0; c < n_channels; ++c) {
			levels->obj[c] = PyTuple_GetItem(levels_obj, c);

			if (levels->all_floats &&
			    PyFloat_Check(levels->obj[c])) {
				const double level_dB =
					PyFloat_AS_DOUBLE(levels->obj[c]);

				levels->fl[c] = dB2lin(level_dB);
			} else {
				levels->all_floats = 0;
			}
		}
	} else {
		PyErr_SetString(PyExc_TypeError,
				"level values must be float or tuple");
		return -1;
	}

	return 0;
}

static int frag_resize(Fragment *frag, size_t length)
{
	const size_t start = frag->length * sizeof(sample_t);
	const size_t size = length * sizeof(sample_t);
	const ssize_t extra = size - start;
	unsigned c;

	for (c = 0; c < frag->n_channels; ++c) {
		if (frag->data[c] == NULL)
			frag->data[c] = PyMem_Malloc(size);
		else
			frag->data[c] = PyMem_Realloc(frag->data[c], size);

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


/* ----------------------------------------------------------------------------
 * _splat methods
 */

PyDoc_STRVAR(splat_lin2dB_doc,
"lin2dB(value)\n"
"\n"
"Convert floating point linear ``value`` to dB.\n");

static PyObject *splat_lin2dB(PyObject *self, PyObject *args)
{
	double level;

	if (!PyArg_ParseTuple(args, "d", &level))
		return NULL;

	return PyFloat_FromDouble(lin2dB(level));
}

PyDoc_STRVAR(splat_dB2lin_doc,
"dB2lin(value)\n"
"\n"
"Convert floating point dB ``value`` to linear.\n");

static PyObject *splat_dB2lin(PyObject *self, PyObject *args)
{
	double dB;

	if (!PyArg_ParseTuple(args, "d", &dB))
		return NULL;

	return PyFloat_FromDouble(dB2lin(dB));
}

/* ----------------------------------------------------------------------------
 * Sources
 */

/* -- helpers for sources -- */

#define splat_check_all_floats(...)			\
	_splat_check_all_floats(PP_NARG(__VA_ARGS__), __VA_ARGS__)

static int _splat_check_all_floats(size_t n, ...)
{
	va_list objs;
	size_t i;

	va_start(objs, n);

	for (i = 0; i < n; ++i) {
		PyObject *obj = va_arg(objs, PyObject *);

		if (!PyFloat_Check(obj))
			break;
	}

	va_end(objs);

	return (i == n);
}

/* -- sine source -- */

static void splat_sine_floats(Fragment *frag, const double *levels,
			      double freq, double phase)
{
	const double k = 2 * M_PI * freq / frag->rate;
	size_t i;

	for (i = 0; i < frag->length; ++i) {
		const double s = sin(k * i);
		unsigned c;

		for (c = 0; c < frag->n_channels; ++c)
			frag->data[c][i] = s * levels[c];
	}
}

static int splat_sine_signals(Fragment *frag, PyObject **levels,
			      PyObject *freq, PyObject *phase)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_AMP,
	};
	static const double k = 2 * M_PI;
	struct splat_signal sig;
	PyObject *signals[SIG_AMP + MAX_CHANNELS];
	unsigned c;

	signals[SIG_FREQ] = freq;
	signals[SIG_PHASE] = phase;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag, signals,
			      (SIG_AMP + frag->n_channels)))
		return -1;

	while (splat_signal_next(&sig) == SIGNAL_VECTOR_CONTINUE) {
		size_t i, j;

		for (i = sig.cur, j = 0; i < sig.end; ++i, j++) {
			const double t = (double)i / frag->rate;
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double s = sin(k * f * (t + ph));

			for (c = 0; c < frag->n_channels; ++c) {
				const double a =
					sig.vectors[SIG_AMP + c].data[j];

				frag->data[c][i] = s * dB2lin(a);
			}
		}
	}

	splat_signal_free(&sig);

	return 0;
}

PyDoc_STRVAR(splat_sine_doc,
"sine(fragment, levels, frequency, phase)\n"
"\n"
"Generate a sine wave for the given ``levels``, ``frequency`` and ``phase`` "
"signals over the entire ``fragment``.\n");

static PyObject *splat_sine(PyObject *self, PyObject *args)
{
	Fragment *frag;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *phase;

	struct splat_levels_helper levels;
	int all_floats;

	if (!PyArg_ParseTuple(args, "O!OOO", &splat_FragmentType, &frag,
			      &levels_obj, &freq, &phase))
		return NULL;

	if (frag_get_levels(frag, &levels, levels_obj))
		return NULL;

	all_floats = levels.all_floats && splat_check_all_floats(freq, phase);

	if (all_floats)
		splat_sine_floats(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				  PyFloat_AS_DOUBLE(phase));
	else if (splat_sine_signals(frag, levels.obj, freq, phase))
		return NULL;

	Py_RETURN_NONE;
}

/* -- square source -- */

static void splat_square_floats(Fragment *frag, const double *fl_pos,
				double freq, double phase, double ratio)
{
	const double k = freq / frag->rate;
	const double ph0 = freq * phase;
	double fl_neg[MAX_CHANNELS];
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

static int splat_square_signals(Fragment *frag, PyObject **levels,
				PyObject *freq, PyObject *phase,
				PyObject *ratio)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_RATIO,
		SIG_AMP,
	};
	struct splat_signal sig;
	PyObject *signals[SIG_AMP + MAX_CHANNELS];
	unsigned c;

	signals[SIG_FREQ] = freq;
	signals[SIG_PHASE] = phase;
	signals[SIG_RATIO] = ratio;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag, signals,
			      (SIG_AMP + frag->n_channels)))
		return -1;

	while (splat_signal_next(&sig) == SIGNAL_VECTOR_CONTINUE) {
		size_t i, j;

		for (i = sig.cur, j = 0; i < sig.end; ++i, ++j) {
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double t = (double)i / frag->rate;
			double n_periods;
			const double t_rel = modf((f * (t + ph)), &n_periods);
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

	return 0;
}

PyDoc_STRVAR(splat_square_doc,
"square(fragment, levels, frequency, phase, ratio=0.5)\n"
"\n"
"Generate a square wave with the given ``ratio`` over the entire "
"``fragment``.\n");

static PyObject *splat_square(PyObject *self, PyObject *args)
{
	Fragment *frag;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *phase;
	PyObject *ratio = splat_init_source_ratio;

	struct splat_levels_helper levels;
	int all_floats;

	if (!PyArg_ParseTuple(args, "O!OOO|O", &splat_FragmentType, &frag,
			      &levels_obj, &freq, &phase, &ratio))
		return NULL;

	if (frag_get_levels(frag, &levels, levels_obj))
		return NULL;

	all_floats = levels.all_floats;
	all_floats = all_floats && splat_check_all_floats(freq, phase, ratio);

	if (all_floats)
		splat_square_floats(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				    PyFloat_AS_DOUBLE(phase),
				    PyFloat_AS_DOUBLE(ratio));
	else if (splat_square_signals(frag, levels.obj, freq, phase, ratio))
		return NULL;

	Py_RETURN_NONE;
}

/* -- triangle source -- */

static void splat_triangle_floats(Fragment *frag, const double *lvls,
				  double freq, double phase, double ratio)
{
	const double k = freq / frag->rate;
	const double ph0 = freq * phase;
	double a1[MAX_CHANNELS], b1[MAX_CHANNELS];
	double a2[MAX_CHANNELS], b2[MAX_CHANNELS];
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

static int splat_triangle_signals(Fragment *frag, PyObject **levels,
				  PyObject *freq, PyObject *phase,
				  PyObject *ratio)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_RATIO,
		SIG_AMP,
	};
	struct splat_signal sig;
	PyObject *signals[SIG_AMP + MAX_CHANNELS];
	unsigned c;

	signals[SIG_FREQ] = freq;
	signals[SIG_PHASE] = phase;
	signals[SIG_RATIO] = ratio;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag, signals,
			      (SIG_AMP + frag->n_channels)))
		return -1;

	while (splat_signal_next(&sig) == SIGNAL_VECTOR_CONTINUE) {
		size_t i, j;

		for (i = sig.cur, j = 0; i < sig.end; ++i, ++j) {
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double t = (double)i / frag->rate;
			double n_periods;
			const double t_rel = modf((f * (t + ph)), &n_periods);
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

	return 0;
}

PyDoc_STRVAR(splat_triangle_doc,
"triangle(fragment, levels, frequency, phase, ratio=0.5)\n"
"\n"
"Generate a triangle wave with the given ``ratio`` over the entire "
"``fragment``.\n");

static PyObject *splat_triangle(PyObject *self, PyObject *args)
{
	Fragment *frag;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *phase;
	PyObject *ratio = splat_init_source_ratio;

	struct splat_levels_helper levels;
	int all_floats;

	if (!PyArg_ParseTuple(args, "O!OOO|O", &splat_FragmentType, &frag,
			      &levels_obj, &freq, &phase, &ratio))
		return NULL;

	if (frag_get_levels(frag, &levels, levels_obj))
		return NULL;

	all_floats = levels.all_floats;
	all_floats = all_floats && splat_check_all_floats(freq, phase, ratio);

	if (all_floats)
		splat_triangle_floats(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				      PyFloat_AS_DOUBLE(phase),
				      PyFloat_AS_DOUBLE(ratio));
	else if (splat_triangle_signals(frag, levels.obj, freq, phase, ratio))
		return NULL;

	Py_RETURN_NONE;
}

/* -- overtones source -- */

struct overtone {
	PyObject *ratio;
	double fl_ratio;
	PyObject *phase;
	double fl_phase;
	struct splat_levels_helper levels;
};

static void splat_overtones_float(Fragment *frag, const double *levels,
				  double freq, double phase,
				  struct overtone *overtones, Py_ssize_t n)
{
	const double k = 2 * M_PI * freq;
	const double max_ratio = (frag->rate / freq) / 2;
	struct overtone *ot;
	const struct overtone *ot_end = &overtones[n];
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

static int splat_overtones_mixed(Fragment *frag, PyObject **levels,
				 PyObject *freq, PyObject *phase,
				 struct overtone *overtones, Py_ssize_t n)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE,
		SIG_AMP,
	};
	const double k = 2 * M_PI;
	const double half_rate = frag->rate / 2;
	PyObject *signals[SIG_AMP + MAX_CHANNELS];
	struct splat_signal sig;
	struct overtone *ot;
	const struct overtone *ot_end = &overtones[n];
	unsigned c;

	signals[SIG_FREQ] = freq;
	signals[SIG_PHASE] = phase;

	for (c = 0; c < frag->n_channels; ++c)
		signals[SIG_AMP + c] = levels[c];

	if (splat_signal_init(&sig, frag, signals,
			      (SIG_AMP + frag->n_channels)))
		return -1;

	while (splat_signal_next(&sig) == SIGNAL_VECTOR_CONTINUE) {
		size_t i, j;

		for (i = sig.cur, j = 0; i < sig.end; ++i, ++j) {
			const double f = sig.vectors[SIG_FREQ].data[j];
			const double max_ratio = half_rate / f;
			const double m = k * f;
			const double ph = sig.vectors[SIG_PHASE].data[j];
			const double t = ph + ((double)i / frag->rate);

			for (ot = overtones; ot != ot_end; ++ot) {
				double s;

				if (ot->fl_ratio >= max_ratio)
					continue;

				s = sin(m * ot->fl_ratio * (t + ot->fl_phase));

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

	return 0;
}

static int splat_overtones_signal(Fragment *frag, PyObject **levels,
				  PyObject *freq, PyObject *phase,
				  struct overtone *overtones, Py_ssize_t n)
{
	enum {
		SIG_FREQ = 0,
		SIG_PHASE = 1,
		SIG_AMP = 2,
		SIG_OT = 2 + MAX_CHANNELS,
	};
	const double k = 2 * M_PI;
	const double half_rate = frag->rate / 2;
	PyObject **signals;
	struct splat_signal sig;
	struct overtone *ot;
	const struct overtone *ot_end = &overtones[n];
	static const size_t sig_freq = 0;
	static const size_t sig_phase = 1;
	static const size_t sig_amp = 2;
	const size_t sig_ot = sig_amp + frag->n_channels;
	/* for each overtone: ratio, phase and levels */
	const size_t sig_n = sig_ot + (n * (2 + frag->n_channels));
	unsigned c;

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

	if (splat_signal_init(&sig, frag, signals, sig_n))
		return -1;

	while (splat_signal_next(&sig) == SIGNAL_VECTOR_CONTINUE) {
		size_t i, j;

		for (i = sig.cur, j = 0; i < sig.end; ++i, ++j) {
			const double t = (double)i / frag->rate;
			const double f = sig.vectors[sig_freq].data[j];
			const double max_ratio = half_rate / f;
			const double ph = sig.vectors[sig_phase].data[j];
			const double m = k * f;
			const double tph = t + ph;
			const struct signal_vector *otv = &sig.vectors[sig_ot];

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

	return 0;
}

PyDoc_STRVAR(splat_overtones_doc,
"overtones(fragment, levels, frequency, phase, overtones)\n"
"\n"
"Generate a sum of overtones as pure sine waves with the given fundamental "
"``frequency`` and ``levels`` in dB.\n"
"\n"
"The ``overtones`` are described with a list of 3-tuples containing the "
"ratio between the overtone and the fundamental frequency, the phase and "
"levels: (ratio, phase, levels).  All these values can be signals, and the "
"levels can either be a single value for all channels or individual values."
"The generation is performed over the entire fragment.\n");

static PyObject *splat_overtones(PyObject *self, PyObject *args)
{
	enum {
		OT_RATIO = 0,
		OT_PHASE,
		OT_LEVELS,
	};
	Fragment *frag;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *phase;
	PyObject *overtones_obj;

	struct splat_levels_helper levels;
	struct overtone *overtones;
	struct overtone *ot;
	Py_ssize_t n;
	Py_ssize_t pos;
	int all_floats;
	int ot_all_floats;
	int stat = 0;

	if (!PyArg_ParseTuple(args, "O!OOOO!", &splat_FragmentType, &frag,
			      &levels_obj, &freq, &phase,
			      &PyList_Type, &overtones_obj))
		return NULL;

	if (frag_get_levels(frag, &levels, levels_obj))
		return NULL;

	n = PyList_GET_SIZE(overtones_obj);
	overtones = PyMem_Malloc(n * sizeof(struct overtone));

	if (overtones == NULL)
		return PyErr_NoMemory();

	all_floats = levels.all_floats && splat_check_all_floats(freq, phase);
	ot_all_floats = 1;

	for (pos = 0, ot = overtones; pos < n; ++pos, ++ot) {
		PyObject *ot_params = PyList_GET_ITEM(overtones_obj, pos);
		PyObject *ot_levels;

		if (!PyTuple_Check(ot_params) ||
		    (PyTuple_GET_SIZE(ot_params) != 3)) {
			PyErr_SetString(PyExc_ValueError,
					"overtone params must be a 3-tuple");
			goto free_overtones;
		}

		ot->ratio = PyTuple_GET_ITEM(ot_params, OT_RATIO);

		if (ot_all_floats && PyFloat_Check(ot->ratio))
			ot->fl_ratio = PyFloat_AS_DOUBLE(ot->ratio);
		else
			ot_all_floats = 0;

		ot->phase = PyTuple_GET_ITEM(ot_params, OT_PHASE);

		if (ot_all_floats && PyFloat_Check(ot->phase))
			ot->fl_phase = PyFloat_AS_DOUBLE(ot->phase);
		else
			ot_all_floats = 0;

		ot_levels = PyTuple_GET_ITEM(ot_params, OT_LEVELS);

		if (frag_get_levels(frag, &ot->levels, ot_levels))
			goto free_overtones;

		ot_all_floats = ot_all_floats && ot->levels.all_floats;
	}

	all_floats = all_floats && ot_all_floats;

	if (all_floats)
		splat_overtones_float(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				      PyFloat_AS_DOUBLE(phase), overtones, n);
	else if (ot_all_floats)
		stat = splat_overtones_mixed(frag, levels.obj, freq, phase,
					     overtones, n);
	else
		stat = splat_overtones_signal(frag, levels.obj, freq, phase,
					      overtones, n);
free_overtones:
	PyMem_Free(overtones);

	if (stat)
		return NULL;

	Py_RETURN_NONE;
}

/* ----------------------------------------------------------------------------
 * Filters
 */

PyDoc_STRVAR(splat_dec_envelope_doc,
"dec_envelope(fragment, k=1.0, p=1.0)\n"
"\n"
"This filter applies a decreasing envelope over the ``fragment`` with ``k`` "
"and ``p`` arguments as follows, for a sound signal ``s`` at index ``i``:\n"
"\n"
".. math::\n"
"\n"
"   s[i] = \\frac{s[i]}{(1 + \\frac{i}{k})^p}\n"
"\n");

static PyObject *splat_dec_envelope(PyObject *self, PyObject *args)
{
	Fragment *frag;
	double k = 1.0;
	double p = 1.0;

	unsigned c;

	if (!PyArg_ParseTuple(args, "O!|dd", &splat_FragmentType, &frag,
			      &k, &p))
		return NULL;

	if (k == 0.0) {
		PyErr_SetString(PyExc_ValueError, "k must not be 0");
		return NULL;
	}

	for (c = 0; c < frag->n_channels; ++c) {
		size_t i;

		for (i = 0; i < frag->length; ++i) {
			const double m = pow(1.0 + ((double)i / k), p);
			frag->data[c][i] /= m;
		}
	}

	Py_RETURN_NONE;
}

PyDoc_STRVAR(splat_reverse_doc,
"reverse(fragment)\n"
"\n"
"Reverse the order of all the ``fragment`` samples.\n");

static PyObject *splat_reverse(PyObject *self, PyObject *args)
{
	Fragment *frag;

	unsigned c;

	if (!PyArg_ParseTuple(args, "O!", &splat_FragmentType, &frag))
		return NULL;

	for (c = 0; c < frag->n_channels; ++c) {
		size_t i;
		size_t j;

		for (i = 0, j = (frag->length - 1); i < j; ++i, --j) {
			const double s = frag->data[c][i];

			frag->data[c][i] = frag->data[c][j];
			frag->data[c][j] = s;
		}
	}

	Py_RETURN_NONE;
}

PyDoc_STRVAR(splat_reverb_doc,
"reverb(fragment, delays, time_factor=0.2, gain_factor=6.0, seed=0)\n"
"\n"
"This filter creates a fast basic reverb effect with some randomness.\n"
"\n"
"The ``delays`` are a list of 2-tuples with a delay duration in seconds and "
"a gain in dB.  They are used to repeat and mix the whole ``fragment`` once "
"for each element in the list, shifted by the given time and amplified by the "
"given gain.  All values must be floating point numbers.  The time delay must "
"not be negative - it's a *causal* reverb.\n"
"\n"
"The ``time_factor`` and ``gain_factor`` parameters are used when adding a "
"random element to the delay and gain.  For example, a ``time_factor`` of 0.2 "
"means the delay will be randomly picked between 1.0 and 1.2 times the value "
"given in the ``delays`` list.  Similarly, for a ``gain_factor`` of 6.0dB the "
"gain will be randomly picked within +/- 6dB around the given value in "
"``delays``.\n"
"\n"
"The ``seed`` argument can be used to initialise the pseudo-random number "
"sequence.  With the default value of 0, the seed will be initialised based "
"on the current time.\n"
"\n"
".. note::\n"
"\n"
"   This filter function can also produce a *delay* effect by specifiying "
"   only a few regularly spaced ``delays``.\n");

static PyObject *splat_reverb(PyObject *self, PyObject *args)
{
	struct delay {
		size_t time;
		double gain;
#if USE_V4SF
		v4sf gain4;
#endif
	};

	Fragment *frag;
	PyObject *delays_list;
	double time_factor = 0.2;
	double gain_factor = 6.0;
	unsigned int seed = 0;

	struct delay *delays[MAX_CHANNELS];
	Py_ssize_t n_delays;
	size_t max_delay;
	size_t max_index;
	size_t d;
	unsigned c;
	size_t i;

	if (!PyArg_ParseTuple(args, "O!O!|ddI", &splat_FragmentType, &frag,
			      &PyList_Type, &delays_list, &time_factor,
			      &gain_factor, &seed))
		return NULL;

	if (!seed)
		seed = time(0);

	srand(seed);

	n_delays = PyList_GET_SIZE(delays_list);

	for (c = 0; c < frag->n_channels; ++c) {
		delays[c] = PyMem_Malloc(n_delays * sizeof(struct delay));

		if (delays[c] == NULL) {
			while (--c)
				PyMem_Free(delays[c]);

			return PyErr_NoMemory();
		}
	}

	max_delay = 0;

	for (d = 0; d < n_delays; ++d) {
		PyObject *pair = PyList_GetItem(delays_list, d);
		double time;
		double gain;

		if (!PyTuple_Check(pair)) {
			PyErr_SetString(PyExc_TypeError,
					"delay values must be a tuple");
			return NULL;
		}

		if (PyTuple_GET_SIZE(pair) != 2) {
			PyErr_SetString(PyExc_ValueError,
					"delay tuple length must be 2");
			return NULL;
		}

		time = PyFloat_AsDouble(PyTuple_GetItem(pair, 0));

		if (time < 0.0) {
			PyErr_SetString(PyExc_ValueError,
					"delay time must be >= 0");
			return NULL;
		}

		for (c = 0; c < frag->n_channels; ++c) {
			const double c_time =
				(time
				 * (1.0 + (rand() * time_factor / RAND_MAX)));
#if USE_V4SF
			delays[c][d].time = c_time * frag->rate / 4;
#else
			delays[c][d].time = c_time * frag->rate;
#endif

			if (delays[c][d].time > max_delay)
				max_delay = delays[c][d].time;
		}

		gain = PyFloat_AsDouble(PyTuple_GetItem(pair, 1));

		for (c = 0; c < frag->n_channels; ++c) {
			const double c_gain_dB =
				(gain - gain_factor
				 + (rand() * gain_factor * 2.0 / RAND_MAX));
#if USE_V4SF
			const double c_gain = dB2lin(c_gain_dB);
			const v4sf gain4 = {c_gain, c_gain, c_gain, c_gain};
			delays[c][d].gain4 = gain4;
			delays[c][d].gain = c_gain;
#else
			delays[c][d].gain = dB2lin(c_gain_dB);
#endif
		}
	}

	max_index = frag->length - 1;

#if USE_V4SF
	max_delay *= 4;
#endif

	if (frag_grow(frag, (frag->length + max_delay)))
		return NULL;

	for (c = 0; c < frag->n_channels; ++c) {
		const struct delay *c_delay = delays[c];
#if USE_V4SF
		v4sf *c_data = (v4sf *)frag->data[c];
#else
		sample_t *c_data = frag->data[c];
#endif

		i = max_index;

#if USE_V4SF
		while (i % 4) {
			const double s = frag->data[c][i];

			for (d = 0; d < n_delays; ++d) {
				const double z = s * c_delay[d].gain;
				frag->data[c][i + (c_delay[d].time * 4)] += z;
			}

			i--;
		}

		i /= 4;
#endif

		do {
#if USE_V4SF
			const v4sf s = c_data[i];
#else
			const double s = c_data[i];
#endif

			for (d = 0; d < n_delays; ++d) {
#if USE_V4SF
				const v4sf z = s * c_delay[d].gain4;
#else
				const double z = s * c_delay[d].gain;
#endif
				c_data[i + c_delay[d].time] += z;
			}
		} while (i--);
	}

	for (c = 0; c < frag->n_channels; ++c)
		PyMem_Free(delays[c]);

	Py_RETURN_NONE;
}

PyDoc_STRVAR(splat_sample_precision_doc,
"sample_precision()\n"
"\n"
"Get the sound sample precision as the number of bits in each floating point "
"sample value.  This is typically 64 or 32 when built with fast mode turned "
"on.\n");

static PyObject *splat_sample_precision(PyObject *self, PyObject *_)
{
	return PyLong_FromLong(sizeof(sample_t) * 8);
}

static PyMethodDef splat_methods[] = {
	{ "lin2dB", splat_lin2dB, METH_VARARGS,
	  splat_lin2dB_doc },
	{ "dB2lin", splat_dB2lin, METH_VARARGS,
	  splat_dB2lin_doc },
	{ "sine", splat_sine, METH_VARARGS,
	  splat_sine_doc },
	{ "square", splat_square, METH_VARARGS,
	  splat_square_doc },
	{ "triangle", splat_triangle, METH_VARARGS,
	  splat_triangle_doc },
	{ "overtones", splat_overtones, METH_VARARGS,
	  splat_overtones_doc },
	{ "dec_envelope", splat_dec_envelope, METH_VARARGS,
	  splat_dec_envelope_doc },
	{ "reverse", splat_reverse, METH_VARARGS,
	  splat_reverse_doc },
	{ "reverb", splat_reverb, METH_VARARGS,
	  splat_reverb_doc },
	{ "sample_precision", splat_sample_precision, METH_NOARGS,
	  splat_sample_precision_doc },
	{ NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC init_splat(void)
{
	struct splat_type {
		PyTypeObject *type;
		const char *name;
	};
	static const struct splat_type splat_types[] = {
		{ &splat_SignalType, "Signal" },
		{ &splat_FragmentType, "Fragment" },
		{ NULL, NULL }
	};
	const struct splat_type *it;
	PyObject *m;

	for (it = splat_types; it->type != NULL; ++it) {
		if (it->type->tp_new == NULL)
			it->type->tp_new = PyType_GenericNew;

		if (PyType_Ready(it->type))
			return;
	}

	m = Py_InitModule("_splat", splat_methods);

	for (it = splat_types; it->type != NULL; ++it) {
		Py_INCREF((PyObject *)it->type);
		PyModule_AddObject(m, it->name, (PyObject *)it->type);
	}

	splat_init_source_ratio = PyFloat_FromDouble(0.5);
	PyModule_AddObject(m, "_init_source_ratio", splat_init_source_ratio);
}

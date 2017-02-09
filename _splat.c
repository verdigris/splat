/*
    Splat - _splat.c

    Copyright (C) 2012, 2013, 2014, 2015, 2016, 2017
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

#define BASE_TYPE_FLAGS (Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE)

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

static const char SPLAT_INT_8[] = "int8";
static const char SPLAT_INT_16[] = "int16";
static const char SPLAT_INT_24[] = "int24";
static const char SPLAT_FLOAT_32[] = "float32";
static const char SPLAT_FLOAT_64[] = "float64";

#define SPLAT_NATIVE_SAMPLE_WIDTH (sizeof(sample_t) * 8)

struct splat_raw_io {
	const char *sample_type;
	size_t sample_width;
	void (*import)(sample_t *out, const char *in, size_t n, size_t step);
	void (*export)(char *out, const sample_t **it, unsigned channels,
		       size_t n);
};

static const struct splat_raw_io *splat_get_raw_io(const char *sample_type);

/* ----------------------------------------------------------------------------
 * Module constants
 */

/* Initial ratio for sound sources (square, triangle) */
static PyObject *splat_init_source_ratio;

/* A dictionary with sample type names as keys and size as values */
static PyObject *splat_sample_types;

/* A Python float with the value of 0.0 */
static PyObject *splat_zero;

/* A Python float with the value of 1.0 */
static PyObject *splat_one;

/* Page utilities */
size_t splat_page_size;
const void *splat_zero_page;

#ifdef SPLAT_FAST
sf_float_t splat_sine_step;
const sf_float_t splat_fast_inc = { 0.0, 1.0, 2.0, 3.0 };
#endif

#if defined(SPLAT_SSE)
const __m128 splat_sse_zero = SPLAT_QUAD(0.0);
const __m128 splat_sse_one = SPLAT_QUAD(1.0);
const __m128 splat_sse_two = SPLAT_QUAD(2.0);
const __m128 splat_sse_pi = SPLAT_QUAD(M_PI);
const __m128 splat_sse_pi2 = SPLAT_QUAD(M_PI * 2.0);
const __m128 splat_sse_inc = { 0.0, 1.0, 2.0, 3.0 };
#endif

int splat_obj2double(PyObject *obj, double *out)
{
	double value;

	if (PyFloat_Check(obj))
		value = PyFloat_AsDouble(obj);
	else if (PyLong_Check(obj))
		value = PyLong_AsDouble(obj);
	else if (PyInt_Check(obj))
		value = (double)PyInt_AsLong(obj);
	else
		return -1;

	if (out != NULL)
		*out = value;

	return 0;
}

/* Levels */

static void splat_levels_init_float(const struct splat_fragment *frag,
				    struct splat_levels *levels,
				    PyObject *levels_obj, double gain)
{
	unsigned c;

	levels->n = frag->n_channels;
	levels->all_floats = 1;

	for (c = 0; c < frag->n_channels; ++c) {
		levels->obj[c] = levels_obj;
		levels->fl[c] = gain;
#ifdef SPLAT_FAST
		levels->flq[c] = sf_set(gain);
#endif
	}
}

static int splat_levels_init_tuple(const struct splat_fragment *frag,
				   struct splat_levels *levels,
				   PyObject *levels_obj)
{
	const Py_ssize_t n_channels = PyTuple_GET_SIZE(levels_obj);
	unsigned c;

	if (n_channels > SPLAT_MAX_CHANNELS) {
		PyErr_SetString(PyExc_ValueError, "too many channels");
		return -1;
	}

	if (n_channels != frag->n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return -1;
	}

	levels->n = n_channels;
	levels->all_floats = 1;

	for (c = 0; c < n_channels; ++c) {
		levels->obj[c] = PyTuple_GetItem(levels_obj, c);

		if (levels->all_floats) {
			double gain;

			if (!splat_obj2double(levels->obj[c], &gain))
				levels->fl[c] = gain;
			else
				levels->all_floats = 0;
		}
	}

	return 0;
}

static void splat_levels_init_signal(const struct splat_fragment *frag,
				     struct splat_levels *levels,
				     PyObject *levels_obj)
{
	unsigned c;

	levels->n = frag->n_channels;
	levels->all_floats = 0;

	for (c = 0; c < frag->n_channels; ++c)
		levels->obj[c] = levels_obj;
}

static int splat_levels_init(const struct splat_fragment *frag,
			     struct splat_levels *levels, PyObject *levels_obj)
{
	double gain;
	int res = 0;

	if (!splat_obj2double(levels_obj, &gain))
		splat_levels_init_float(frag, levels, levels_obj, gain);
	else if (PyTuple_Check(levels_obj))
		res = splat_levels_init_tuple(frag, levels, levels_obj);
	else
		splat_levels_init_signal(frag, levels, levels_obj);

	return res;
}

/* Fragment class declaration */

struct Fragment_object {
	PyObject_HEAD;
	int init;
	struct splat_fragment frag;
};
typedef struct Fragment_object Fragment;

static PyTypeObject splat_FragmentType;

/* ----------------------------------------------------------------------------
 * Spline class
 */

struct Spline_object {
	PyObject_HEAD;
	int init;
	struct splat_spline spline;
};
typedef struct Spline_object Spline;

static PyTypeObject splat_SplineType;

struct splat_spline *splat_spline_from_obj(PyObject *obj)
{
	if (!PyObject_TypeCheck(obj, &splat_SplineType))
		return NULL;

	return &((Spline *)obj)->spline;
}

static void Spline_dealloc(Spline *self)
{
	if (self->init) {
		Py_DECREF(self->spline.pols);
		self->init = 0;
	}
}

static int Spline_init(Spline *self, PyObject *args)
{
	PyObject *poly;
	struct splat_spline *spline = &self->spline;
	PyObject *db;

	if (!PyArg_ParseTuple(args, "O!dO!", &PyList_Type, &spline->pols,
			      &spline->k0, &PyBool_Type, &db))
		return -1;

	Py_INCREF(spline->pols);

	poly = PyList_GET_ITEM(spline->pols, 0);
	spline->start = PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(poly, 0));

	poly = PyList_GET_ITEM(spline->pols, PyList_GET_SIZE(spline->pols)-1);
	spline->end = PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(poly, 1));

	spline->db = (db == Py_True) ? 1 : 0;

	self->init = 1;

	return 0;
}

static PyObject *Spline_new(PyTypeObject *type, PyObject *args, PyObject *kw)
{
	Spline *self;

	self = (Spline *)type->tp_alloc(type, 0);

	if (self == NULL)
		return PyErr_NoMemory();

	self->init = 0;

	return (PyObject *)self;
}

static PyTypeObject splat_SplineType = {
	PyObject_HEAD_INIT(NULL)
	0,                                 /* ob_size */
	"_splat.Spline",                   /* tp_name */
	sizeof(Spline),                    /* tp_basicsize */
	0,                                 /* tp_itemsize */
	(destructor)Spline_dealloc,        /* tp_dealloc */
	0,                                 /* tp_print */
	0,                                 /* tp_getattr */
	0,                                 /* tp_setattr */
	0,                                 /* tp_compare */
	0,                                 /* tp_repr */
	0,                                 /* tp_as_number */
	0,                                 /* tp_as_sequence */
	0,                                 /* tp_as_mapping */
	0,                                 /* tp_hash  */
	0,                                 /* tp_call */
	0,                                 /* tp_str */
	0,                                 /* tp_getattro */
	0,                                 /* tp_setattro */
	0,                                 /* tp_as_buffer */
	BASE_TYPE_FLAGS,                   /* tp_flags */
	0,                                 /* tp_doc */
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
	(initproc)Spline_init,             /* tp_init */
	0,                                 /* tp_alloc */
	Spline_new,                        /* tp_new */
};

/* -- Signal class -- */

struct Signal_object {
	PyObject_HEAD;
	int init;
	PyObject **signals;
	struct splat_signal sig;
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
		self->init = 0;
	}

	self->ob_type->tp_free((PyObject *)self);
}

static int Signal_init(Signal *self, PyObject *args)
{
	Fragment *frag;
	PyObject *sig_obj;
	double duration = 0.0;
	double origin = 0.0;

	size_t n_signals;
	size_t i;
	size_t length;

	if (!PyArg_ParseTuple(args, "O!O|dd", &splat_FragmentType, &frag,
			      &sig_obj, &duration, &origin))
		return -1;

	if (PyTuple_Check(sig_obj))
		n_signals = PyTuple_GET_SIZE(sig_obj);
	else
		n_signals = 1;

	if (duration != 0.0) {
		if (duration < 0.0) {
			PyErr_SetString(PyExc_ValueError, "negative duration");
			return -1;
		}

		length = duration * frag->frag.rate;
	} else {
		length = frag->frag.length;
	}

	if (origin < 0.0) {
		PyErr_SetString(PyExc_ValueError, "negative signal origin");
		return -1;
	}

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

	if (splat_signal_init(&self->sig, length, (origin * frag->frag.rate),
			      self->signals, n_signals, frag->frag.rate)) {
		Signal_free_signals(self, n_signals);
		return -1;
	}

	self->offset = 0;
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
	return self->sig.length;
}

static PyObject *SignalObj_sq_item(Signal *self, Py_ssize_t i)
{
	ssize_t offset;

	offset = splat_signal_get(&self->sig, i);

	if (offset < 0) {
		PyErr_SetString(PyExc_IndexError, "out of signal range");
		return NULL;
	}

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
	0,                                 /* tp_doc */
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

/* ----------------------------------------------------------------------------
 * Fragment class
 */

static PyObject *splat_frag_peak_as_dict(const struct splat_peak *peak)
{
	return Py_BuildValue("{sdsdsdsd}", "avg", peak->avg, "max", peak->max,
			     "min", peak->min, "peak", peak->peak);
}

struct splat_fragment *splat_frag_from_obj(PyObject *obj)
{
	if (!PyObject_TypeCheck(obj, &splat_FragmentType))
		return NULL;

	return &((Fragment *)obj)->frag;
}

static void Fragment_dealloc(Fragment *self)
{
	if (self->init) {
		splat_frag_free(&self->frag);
		self->init = 0;
	}

	self->ob_type->tp_free((PyObject *)self);
}

static int splat_frag_mmap(struct splat_fragment *frag, unsigned n_channels,
			   unsigned rate, size_t length, const char *name,
			   PyObject *obj)
{
	const char *new_path;

	if (obj == Py_True) {
		new_path = NULL;
	} else if ((obj != NULL) && PyString_Check(obj)) {
		new_path = PyString_AsString(obj);
	} else {
		PyErr_SetString(PyExc_ValueError,
				"invalid mmap argument");
		return -1;
	}

	return splat_frag_init_mmap(frag, n_channels, rate, length, name,
				    new_path);
}

static int Fragment_init(Fragment *self, PyObject *args, PyObject *kw)
{
	static char *kwlist[] = {
		"channels", "rate", "duration", "length", "name", "mmap",
		NULL };
	unsigned n_channels = 2;
	unsigned rate = 48000;
	double duration = 0.0;
	unsigned long length = 0;
	const char *name = NULL;
	PyObject *mmap_obj = NULL;

	int ret;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "|IIdkzO", kwlist,
					 &n_channels, &rate, &duration,
					 &length, &name, &mmap_obj))
		return -1;

	if (n_channels > SPLAT_MAX_CHANNELS) {
		PyErr_SetString(PyExc_ValueError,
				"exceeding maximum number of channels");
		return -1;
	}

	if (!rate) {
		PyErr_SetString(PyExc_ValueError, "rate cannot be 0 Hz");
		return -1;
	}

	if (duration < 0.0) {
		PyErr_SetString(PyExc_ValueError, "negative duration");
		return -1;
	}

	if (!length) {
		length = duration * rate;
	} else if (duration != 0.0) {
		PyErr_SetString(PyExc_ValueError,
				"cannot specify both length and duration");
		return -1;
	}

	if (mmap_obj == NULL)
		ret = splat_frag_init(&self->frag, n_channels, rate, length,
				      name);
	else
		ret = splat_frag_mmap(&self->frag, n_channels, rate, length,
				      name, mmap_obj);

	if (ret)
		return -1;

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
	return self->frag.length;
}

static PyObject *Fragment_sq_item(Fragment *self, Py_ssize_t i)
{
	PyObject *sample;
	unsigned c;

	if ((i < 0) || (i >= self->frag.length)) {
		PyErr_SetString(PyExc_IndexError, "index out of range");
		return NULL;
	}

	sample = PyTuple_New(self->frag.n_channels);

	if (sample == NULL)
		return NULL;

	for (c = 0; c < self->frag.n_channels; ++c) {
		const sample_t s = self->frag.channels[c].data[i];
		PyTuple_SET_ITEM(sample, c, PyFloat_FromDouble(s));
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

	if (PyTuple_GET_SIZE(v) != self->frag.n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return -1;
	}

	if ((i < 0) || (i >= self->frag.length)) {
		PyErr_SetString(PyExc_IndexError, "set index error");
		return -1;
	}

	for (c = 0; c < self->frag.n_channels; ++c) {
		PyObject *s = PyTuple_GET_ITEM(v, c);

		if (!PyFloat_CheckExact(s)) {
			PyErr_SetString(PyExc_TypeError,
					"item must contain floats");
			return -1;
		}

		self->frag.channels[c].data[i] =
			(sample_t)PyFloat_AS_DOUBLE(s);
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

PyDoc_STRVAR(rate_doc, "Get the sample rate in Hz.");

static PyObject *Fragment_get_rate(Fragment *self, void *_)
{
	return Py_BuildValue("I", self->frag.rate);
}

PyDoc_STRVAR(duration_doc, "Get the fragment duration in seconds.");

static PyObject *Fragment_get_duration(Fragment *self, void *_)
{
	return Py_BuildValue("f", (double)self->frag.length / self->frag.rate);
}

PyDoc_STRVAR(channels_doc, "Get the number of channels.");

static PyObject *Fragment_get_channels(Fragment *self, void *_)
{
	return Py_BuildValue("I", self->frag.n_channels);
}

PyDoc_STRVAR(name_doc, "Get and set the fragment name.");

static PyObject *Fragment_get_name(Fragment *self, void *_)
{
	if (self->frag.name == NULL)
		Py_RETURN_NONE;

	return PyString_FromString(self->frag.name);
}

static int Fragment_set_name(Fragment *self, PyObject *value, void *_)
{
	if (!PyString_Check(value)) {
		PyErr_SetString(PyExc_TypeError,
				"Fragment name must be a string");
		return -1;
	}

	return splat_frag_set_name(&self->frag, PyString_AS_STRING(value));
}

static PyGetSetDef Fragment_getsetters[] = {
	{ "rate", (getter)Fragment_get_rate, NULL, rate_doc },
	{ "duration", (getter)Fragment_get_duration, NULL, duration_doc },
	{ "channels", (getter)Fragment_get_channels, NULL, channels_doc },
	{ "name", (getter)Fragment_get_name, (setter)Fragment_set_name,
	  name_doc },
	{ NULL }
};

/* Fragment methods */

static void splat_import_float64(sample_t *out, const char *in, size_t n,
				 size_t step)
{
	while (n--) {
		*out++ = *(const double *)in;
		in += step;
	}
}

static void splat_export_float64(char *out, const sample_t **in,
				 unsigned channels, size_t n)
{
	double *out64 = (double *)out;

	while (n--) {
		unsigned c;

		for (c = 0; c < channels; ++c)
			*out64++ = *(in[c]++);
	}
}

static void splat_import_float32(sample_t *out, const char *in, size_t n,
				 size_t step)
{
	while (n--) {
		*out++ = *(const float *)in;
		in += step;
	}
}

static void splat_export_float32(char *out, const sample_t **in,
				 unsigned channels, size_t n)
{
	float *out32 = (float *)out;

	while (n--) {
		unsigned c;

		for (c = 0; c < channels; ++c)
			*out32++ = *(in[c]++);
	}
}

static void splat_import_int24(sample_t *out, const char *in, size_t n,
			       size_t step)
{
	static const int32_t neg24 = 1 << 23;
	static const int32_t neg_mask32 = 0xFF000000;
	static const sample_t scale = (1 << 23) - 1;

	while (n--) {
		const uint8_t *in8 = (const uint8_t *)in;
		int32_t sample32;

		sample32 = *in8++;
		sample32 += (*in8++) << 8;
		sample32 += (*in8++) << 16;

		if (sample32 & neg24)
			sample32 |= neg_mask32;

		*out++ = sample32 / scale;
		in += step;
	}
}

static void splat_export_int24(char *out, const sample_t **in,
			       unsigned channels, size_t n)
{
	static const int32_t scale = (1 << 23) - 1;

	while (n--) {
		unsigned c;

		for (c = 0; c < channels; ++c) {
			const sample_t z = *(in[c]++);
			int32_t s;

			if (z < -1.0)
				s = -scale;
			else if (z > 1.0)
				s = scale;
			else
				s = z * scale;

			*out++ = s & 0xFF;
			*out++ = (s >> 8) & 0xFF;
			*out++ = (s >> 16) & 0xFF;
		}
	}
}

static void splat_import_int16(sample_t *out, const char *in, size_t n,
			       size_t step)
{
	static const sample_t scale = (1 << 15) - 1;

	while (n--) {
		*out++ = *(int16_t *)in / scale;
		in += step;
	}
}

static void splat_export_int16(char *out, const sample_t **in,
			       unsigned channels, size_t n)
{
	static const long scale = (1 << 15) - 1;

	while (n--) {
		unsigned c;

		for (c = 0; c < channels; ++c) {
			const sample_t z = *(in[c]++);
			int16_t s;

			if (z < -1.0)
				s = -scale;
			else if (z > 1.0)
				s = scale;
			else
				s = z * scale;

			*out++ = s & 0xFF;
			*out++ = (s >> 8) & 0xFF;
		}
	}
}

static void splat_import_int8(sample_t *out, const char *in, size_t n,
			      size_t step)
{
	static const sample_t scale = 127.0;

	while (n--) {
		*out++ = *(int8_t *)in / scale;
		in += step;
	}
}

static void splat_export_int8(char *out, const sample_t **in,
			      unsigned channels, size_t n)
{
	while (n--) {
		unsigned c;

		for (c = 0; c < channels; ++c) {
			const sample_t z = (*in[c]++);
			int8_t s;

			if (z < -1.0)
				s = 0;
			else if (z > 1.0)
				s = 255;
			else
				s = z * 127.0;

			*out++ = s;
		}
	}
}

static const struct splat_raw_io splat_raw_io_table[] = {
	{ SPLAT_FLOAT_64, 64, splat_import_float64, splat_export_float64 },
	{ SPLAT_FLOAT_32, 32, splat_import_float32, splat_export_float32 },
	{ SPLAT_INT_24, 24, splat_import_int24, splat_export_int24, },
	{ SPLAT_INT_16, 16, splat_import_int16, splat_export_int16, },
	{ SPLAT_INT_8, 8, splat_import_int8, splat_export_int8 },
};

static const struct splat_raw_io *splat_get_raw_io(const char *sample_type)
{
	const struct splat_raw_io *it;
	const struct splat_raw_io * const end =
		&splat_raw_io_table[ARRAY_SIZE(splat_raw_io_table)];

	for (it = splat_raw_io_table; it != end; ++it)
		if (!strcmp(sample_type, it->sample_type))
			return it;

	PyErr_SetString(PyExc_ValueError, "unsupported sample type");

	return NULL;
}

PyDoc_STRVAR(Fragment_import_bytes_doc,
"import_bytes(raw_bytes, rate, channels, sample_type=splat.SAMPLE_TYPE, "
"offset=None, start=None, end=None)\n"
"\n"
"Import data as raw bytes.\n"
"\n"
"The ``sample_type`` gives the format of the raw data to import as samples, "
":ref:`sample_formats` for more details. "
"The ``rate`` and ``channels`` need to match the Fragment instance values. "
"The ``offset`` argument can be used as a sample number to specify the point "
"where the data starts to be imported into the fragment. "
"The ``start`` and ``end`` arguments can be defined as sample numbers to "
"import only a subset of the data.\n");

static PyObject *Fragment_import_bytes(Fragment *self, PyObject *args,
				       PyObject *kw)
{
	static char *kwlist[] = {
		"raw_bytes", "rate", "channels", "sample_type",
		"offset", "start", "end", NULL };
	PyObject *bytes_obj;
	unsigned rate;
	unsigned n_channels;
	const char *sample_type = SPLAT_NATIVE_SAMPLE_TYPE;
	PyObject *offset_obj = Py_None;
	PyObject *start_obj = Py_None;
	PyObject *end_obj = Py_None;

	size_t sample_size;
	size_t bytes_size;
	size_t frame_size;
	size_t bytes_length;
	const struct splat_raw_io *io;
	size_t offset;
	size_t start;
	size_t end;
	size_t length;
	const char *bytes;
	unsigned c;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!II|sOOO", kwlist,
					 &PyByteArray_Type, &bytes_obj, &rate,
					 &n_channels, &sample_type,
					 &offset_obj, &start_obj, &end_obj))
		return NULL;

	if (rate != self->frag.rate) {
		PyErr_SetString(PyExc_ValueError, "wrong sample rate");
		return NULL;
	}

	if (n_channels != self->frag.n_channels) {
		PyErr_SetString(PyExc_ValueError, "wrong number of channels");
		return NULL;
	}

	io = splat_get_raw_io(sample_type);

	if (io == NULL)
		return NULL;

	bytes_size = PyByteArray_Size(bytes_obj);
	sample_size = io->sample_width / 8;
	frame_size = sample_size * n_channels;
	bytes_length = bytes_size / frame_size;

	if (bytes_size % frame_size) {
		PyErr_SetString(PyExc_ValueError,
				"buffer length not multiple of frame size");
		return NULL;
	}

	if (offset_obj == Py_None)
		offset = 0;
	else if (splat_frag_sample_number(&offset, 0, LONG_MAX, offset_obj))
		return NULL;

	if (start_obj == Py_None)
		start = 0;
	else if (splat_frag_sample_number(&start, 0, bytes_length, start_obj))
		return NULL;

	if (end_obj == Py_None)
		end = bytes_length;
	else if (splat_frag_sample_number(&end, start, bytes_length, end_obj))
		return NULL;

	length = end - start;

	if (splat_frag_grow(&self->frag, (offset + length)))
		return NULL;

	bytes = PyByteArray_AsString(bytes_obj);

	for (c = 0; c < self->frag.n_channels; ++c) {
		const char *in =
			bytes + (start * frame_size) + (c * sample_size);
		sample_t *out = &self->frag.channels[c].data[offset];

		io->import(out, in, length, frame_size);
	}

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_export_bytes_doc,
"export_bytes(sample_type=splat.SAMPLE_TYPE, start=None, end=None)\n"
"\n"
"Export audio data as raw bytes.\n"
"\n"
"The ``sample_type`` is a string to specify the format of the exported "
"samples.  See :ref:`sample_formats` for more details.\n"
"\n"
"The ``start`` and ``end`` arguments can be specified in sample numbers to "
"only get a subset of the data.\n"
"\n"
"A ``bytearray`` object is then returned with the exported sample data.\n");

static PyObject *Fragment_export_bytes(Fragment *self, PyObject *args,
				       PyObject *kw)
{
	static char *kwlist[] = { "sample_type", "start", "end", NULL };
	const char *sample_type = SPLAT_NATIVE_SAMPLE_TYPE;
	PyObject *start_obj = Py_None;
	PyObject *end_obj = Py_None;

	struct splat_fragment *frag = &self->frag;
	const struct splat_raw_io *io;
	size_t sample_size;
	size_t frame_size;
	size_t start;
	size_t end;
	size_t length;
	unsigned c;
	PyObject *bytes_obj;
	Py_ssize_t bytes_size;
	char *out;
	const sample_t *in[SPLAT_MAX_CHANNELS];

	if (!PyArg_ParseTupleAndKeywords(args, kw, "|sOO", kwlist,
					 &sample_type, &start_obj, &end_obj))
		return NULL;

	io = splat_get_raw_io(sample_type);

	if (io == NULL)
		return NULL;

	sample_size = io->sample_width / 8;
	frame_size = sample_size * frag->n_channels;

	if (start_obj == Py_None)
		start = 0;
	else if (splat_frag_sample_number(&start, 0, frag->length, start_obj))
		return NULL;

	if (end_obj == Py_None)
		end = frag->length;
	else if (splat_frag_sample_number(&end, start, frag->length, end_obj))
		return NULL;

	length = end - start;
	bytes_size = length * frame_size;
	bytes_obj = PyByteArray_FromStringAndSize(NULL, bytes_size);

	if (bytes_obj == NULL)
		return PyErr_NoMemory();

	out = PyByteArray_AS_STRING(bytes_obj);

	for (c = 0; c < frag->n_channels; ++c)
		in[c] = &frag->channels[c].data[start];

	io->export(out, in, frag->n_channels, length);

	return bytes_obj;
}

PyDoc_STRVAR(Fragment_mix_doc,
"mix(fragment, offset=0.0, skip=0.0, levels=None, duration=None)\n"
"\n"
"Mix the given other ``fragment`` data into this instance.\n"
"\n"
"This is achieved by simply adding the corresponding samples of an incoming "
"fragment to this fragment' samples.  The ``offset``, ``skip`` and "
"``duration`` values in seconds can be used to alter the mixing times.  The "
"incoming fragment can start being mixed with an ``offset`` into this "
"fragment, its beginning can be skipped until the given ``skip`` time, and "
"the ``duration`` to be mixed can be manually limited.  These values will be "
"automatically adjusted to remain within the available incoming data.  The "
"length of this fragment will be automatically increased if necessary to hold "
"the mixed data.\n"
"\n"
"The ``levels`` argument can be used to alter the amplitude of the incoming "
"fragment while mixing - this does not affect the original fragment.\n"
"\n"
"Please note that the two fragments must use the same sample rate and have "
"the same number of channels.\n");

static PyObject *Fragment_mix(Fragment *self, PyObject *args, PyObject *kw)
{
	static char *kwlist[] = {
		"frag", "offset", "skip", "levels", "duration", NULL };
	Fragment *incoming_obj;
	double offset = 0.0;
	double skip = 0.0;
	PyObject *levels_obj = Py_None;
	PyObject *duration_obj = Py_None;

	struct splat_fragment *frag = &self->frag;
	const struct splat_fragment *incoming;
	struct splat_levels levels;
	ssize_t length;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "O!|ddOO", kwlist,
					 &splat_FragmentType, &incoming_obj,
					 &offset, &skip, &levels_obj,
					 &duration_obj))
		return NULL;

	incoming = &incoming_obj->frag;

	if (incoming->n_channels != frag->n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return NULL;
	}

	if (incoming->rate != frag->rate) {
		PyErr_SetString(PyExc_ValueError, "sample rate mismatch");
		return NULL;
	}

	if (levels_obj == Py_None)
		levels_obj = splat_zero;

	if (splat_levels_init(frag, &levels, levels_obj))
		return NULL;

	if (duration_obj != Py_None) {
		if (!PyFloat_Check(duration_obj)) {
			PyErr_SetString(PyExc_ValueError,
					"duration must be float");
			return NULL;
		}

		length = PyFloat_AS_DOUBLE(duration_obj) * frag->rate;
	} else {
		length = incoming->length;
	}

	if (splat_frag_mix(frag, incoming, &levels, length, offset, skip,
			   levels_obj == splat_zero))
		return NULL;

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_get_peak_doc,
"get_peak()\n"
"\n"
"Get peak and average values for the whole fragment and for each channel.\n"
"\n"
"Scan all the data and look for the peak maximum, minimum and absolute values "
"as well as the average values for each channel and for the whole fragment. "
"The results are returned as a 2-tuple, the first item being for the whole "
"fragment and the second one a list with results for each channel.  "
"Each result is a dictionary with ``max``, ``min``, ``peak`` and ``avg`` "
"values respectively for linear maximum, minimum, absolute peak and average "
"values.\n");

static PyObject *Fragment_get_peak(Fragment *self, PyObject *_)
{
	struct splat_peak chan_peak[SPLAT_MAX_CHANNELS];
	struct splat_peak frag_peak;
	unsigned c;

	PyObject *ret;
	PyObject *frag_peak_obj;
	PyObject *chan_peak_obj;

	splat_frag_get_peak(&self->frag, chan_peak, &frag_peak, 1);
	frag_peak_obj = splat_frag_peak_as_dict(&frag_peak);

	if (frag_peak_obj == NULL)
		return PyErr_NoMemory();

	chan_peak_obj = PyList_New(self->frag.n_channels);

	if (chan_peak_obj == NULL)
		return PyErr_NoMemory();

	for (c = 0; c < self->frag.n_channels; ++c) {
		PyObject *peak_dict = splat_frag_peak_as_dict(&chan_peak[c]);

		if (peak_dict == NULL)
			return PyErr_NoMemory();

		PyList_SET_ITEM(chan_peak_obj, c, peak_dict);
	}

	ret = Py_BuildValue("(OO)", frag_peak_obj, chan_peak_obj);
	Py_DECREF(frag_peak_obj);
	Py_DECREF(chan_peak_obj);

	return ret;
}

PyDoc_STRVAR(Fragment_normalize_doc,
"normalize(level=-0.05, zero=True)\n"
"\n"
"Normalize the amplitude.\n"
"\n"
"The ``level`` value in dB is the resulting maximum amplitude after "
"normalization.  The default value of -0.05 dB is the maximum level while "
"ensuring no clipping occurs due to rounding errors, even when saving with "
"8-bit sample resolution. "
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
	double level_dB = -0.05;
	PyObject *zero = NULL;

	int do_zero;

	if (!PyArg_ParseTuple(args, "|dO!", &level_dB, &PyBool_Type, &zero))
		return NULL;

	if (self->frag.n_channels > SPLAT_MAX_CHANNELS) {
		PyErr_SetString(PyExc_ValueError, "too many channels");
		return NULL;
	}

	do_zero = ((zero == NULL) || (zero == Py_True)) ? 1 : 0;
	splat_frag_normalize(&self->frag, level_dB, do_zero);

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_amp_doc,
"amp(gain)\n"
"\n"
"Amplify the fragment by the given linear ``gain`` which can either be a "
"single signal  to apply to all channels or a tuple with a signal "
"for each individual channel.  Please note that a negative value will "
"invert the data in the fragment.\n");

static PyObject *Fragment_amp(Fragment *self, PyObject *args)
{
	PyObject *gain_obj;

	struct splat_levels gains;

	if (!PyArg_ParseTuple(args, "O", &gain_obj))
		return NULL;

	if (splat_levels_init(&self->frag, &gains, gain_obj))
		return NULL;

	if (splat_frag_amp(&self->frag, &gains))
		return NULL;

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_lin2dB_doc,
"lin2dB()\n"
"\n"
"Convert all the sample values to dB.\n");

static PyObject *Fragment_lin2dB(Fragment *self, PyObject *_)
{
	splat_frag_lin2dB(&self->frag);

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_dB2lin_doc,
"dB2lin()\n"
"\n"
"Convert all the sample values from dB to linear scale.\n");

static PyObject *Fragment_dB2lin(Fragment *self, PyObject *_)
{
	splat_frag_dB2lin(&self->frag);

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_offset_doc,
"offset(value, start=0.0)\n"
"\n"
"Add a linear offset to the data already in the fragment starting at the "
"``start`` time in seconds.  This is especially useful when generating a "
"modulation fragment.\n");

static PyObject *Fragment_offset(Fragment *self, PyObject *args)
{
	PyObject *offset;
	double start = 0.0;

	if (!PyArg_ParseTuple(args, "O|d", &offset, &start))
		return NULL;

	if (splat_frag_offset(&self->frag, offset, start))
		return NULL;

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_resize_doc,
"resize(duration=0.0, length=0)\n"
"\n"
"Resize the fragment to the given ``duration`` in seconds or ``length`` in "
"number of samples per channel.  If the fragment grows, silence is added at "
"the end.  When shrinking, the end of the fragment is lost.\n");

static PyObject *Fragment_resize(Fragment *self, PyObject *args, PyObject *kw)
{
	static char *kwlist[] = { "duration", "length", NULL };
	double duration = 0.0;
	unsigned long length = 0;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "|dk", kwlist,
					 &duration, &length))
		return NULL;

	if (duration < 0.0) {
		PyErr_SetString(PyExc_ValueError, "negative duration");
		return NULL;
	}

	if (!length) {
		length = duration * self->frag.rate;
	} else if (duration != 0.0) {
		PyErr_SetString(PyExc_ValueError,
				"cannot specify both length and duration");
		return NULL;
	}

	if (splat_frag_resize(&self->frag, length))
		return NULL;

	Py_RETURN_NONE;
}

PyDoc_STRVAR(Fragment_resample_doc,
"resample(rate=None, ratio=1.0)\n"
"\n"
"Resample the fragment with a new sample ``rate`` and stretch the "
"fragment duration by the given ``ratio``.  When no ``rate`` is "
"provided, the current sample rate is preserved.  By default, "
"``ratio`` equals to 1.0 so the fragment duration remains unchanged "
"when only resampling with a different rate.  While ``rate`` needs to "
"be an integer as a fragment's sample rate is fixed, ``ratio`` can be "
"a signal to create modulations and effects.\n"
"\n"
"It is also worth noting that no filter is being applied by "
"this function, so down-sampling to a lower rate or using a ``ratio`` "
"smaller than 1.0 may cause spectrum overlap if the input signal contains "
"frequencies higher than half of the new rate (Nyquist theorem).\n");

static PyObject *Fragment_resample(Fragment *self, PyObject *args, PyObject *kw)
{
	struct splat_fragment *frag = &self->frag;

	static char *kwlist[] = { "rate", "ratio", NULL };
	unsigned rate = frag->rate;
	PyObject *ratio = splat_one;

	struct splat_fragment old_frag;
	unsigned c;
	int res;

	if (!PyArg_ParseTupleAndKeywords(args, kw, "|IO", kwlist,
					 &rate, &ratio))
		return NULL;

	if (!rate) {
		PyErr_SetString(PyExc_ValueError, "rate cannot be 0 Hz");
		return NULL;
	}

	if (frag->length < 3) {
		PyErr_SetString(PyExc_ValueError, "fragment is too short");
		return NULL;
	}

	if (splat_frag_init(&old_frag, frag->n_channels, frag->rate,
			    frag->length, NULL))
		return NULL;

	for (c = 0; c < frag->n_channels; ++c) {
		sample_t *mv;

		mv = old_frag.channels[c].data;
		old_frag.channels[c].data = frag->channels[c].data;
		frag->channels[c].data = mv;
	}

	res = splat_frag_resample(frag, &old_frag, rate, ratio);
	splat_frag_free(&old_frag);

	if (res)
		return NULL;

	Py_RETURN_NONE;
}

static PyMethodDef Fragment_methods[] = {
	{ "import_bytes", (PyCFunction)Fragment_import_bytes, METH_KEYWORDS,
	  Fragment_import_bytes_doc },
	{ "export_bytes", (PyCFunction)Fragment_export_bytes, METH_KEYWORDS,
	  Fragment_export_bytes_doc },
	{ "mix", (PyCFunction)Fragment_mix, METH_KEYWORDS,
	  Fragment_mix_doc },
	{ "get_peak", (PyCFunction)Fragment_get_peak, METH_NOARGS,
	  Fragment_get_peak_doc },
	{ "normalize", (PyCFunction)Fragment_normalize, METH_VARARGS,
	  Fragment_normalize_doc },
	{ "amp", (PyCFunction)Fragment_amp, METH_VARARGS,
	  Fragment_amp_doc },
	{ "lin2dB", (PyCFunction)Fragment_lin2dB, METH_NOARGS,
	  Fragment_lin2dB_doc },
	{ "dB2lin", (PyCFunction)Fragment_dB2lin, METH_NOARGS,
	  Fragment_dB2lin_doc },
	{ "offset", (PyCFunction)Fragment_offset, METH_VARARGS,
	  Fragment_offset_doc },
	{ "resize", (PyCFunction)Fragment_resize, METH_KEYWORDS,
	  Fragment_resize_doc },
	{ "resample", (PyCFunction)Fragment_resample, METH_KEYWORDS,
	  Fragment_resample_doc },
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

PyDoc_STRVAR(splat_gen_ref_doc,
"gen_ref(frag)\n"
"\n"
"Generate a reference signal into a single channel fragment."
"\n"
"This is useful mainly for benchmarks and testing purposes.\n");

static PyObject *splat_gen_ref(PyObject *self, PyObject *args)
{
	Fragment *frag_obj;
	struct splat_fragment *frag;
	sample_t *data;
	size_t n;

	if (!PyArg_ParseTuple(args, "O!", &splat_FragmentType, &frag_obj))
		return NULL;

	frag = &frag_obj->frag;

	if (frag->n_channels != 1) {
		PyErr_SetString(PyExc_ValueError,
				"fragment must have a single channel");
		return NULL;
	}

	if (frag->length == 0) {
		PyErr_SetString(PyExc_ValueError,
				"fragment length must be greater than 0");
		return NULL;
	}

	data = frag->channels[0].data;
	n = frag->length;

	while (n--)
		*data++ = n;

	Py_RETURN_NONE;
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

PyDoc_STRVAR(splat_sine_doc,
"sine(fragment, levels, frequency, phase=0.0, origin=0.0)\n"
"\n"
"Generate a sine wave for the given ``levels``, ``frequency`` and ``phase`` "
"signals over the entire ``fragment`` with the given ``origin`` in time.\n");

static PyObject *splat_sine(PyObject *self, PyObject *args)
{
	Fragment *frag_obj;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *phase = splat_zero;
	double origin = 0.0;

	struct splat_fragment *frag;
	struct splat_levels levels;
	int all_floats;

	if (!PyArg_ParseTuple(args, "O!OO|Od", &splat_FragmentType, &frag_obj,
			      &levels_obj, &freq, &phase, &origin))
		return NULL;

	frag = &frag_obj->frag;

	if (splat_levels_init(frag, &levels, levels_obj))
		return NULL;

	all_floats = levels.all_floats && splat_check_all_floats(freq, phase);

	if (all_floats) {
		origin = fmod(origin, 1 / PyFloat_AS_DOUBLE(freq));
		splat_sine_floats(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				  PyFloat_AS_DOUBLE(phase) + origin);
	} else if (splat_sine_signals(frag, levels.obj, freq, phase, origin)) {
		return NULL;
	}

	Py_RETURN_NONE;
}

PyDoc_STRVAR(splat_square_doc,
"square(fragment, levels, frequency, phase=0.0, origin=0.0, ratio=0.5)\n"
"\n"
"Generate a square wave with standard argments and the given ``ratio`` over "
"the entire ``fragment``.  The ratio is between the duration of the high and "
"low states.\n");

static PyObject *splat_square(PyObject *self, PyObject *args)
{
	Fragment *frag_obj;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *phase = splat_zero;
	double origin = 0.0;
	PyObject *ratio = splat_init_source_ratio;

	struct splat_fragment *frag;
	struct splat_levels levels;
	int all_floats;

	if (!PyArg_ParseTuple(args, "O!OO|OdO", &splat_FragmentType, &frag_obj,
			      &levels_obj, &freq, &phase, &origin, &ratio))
		return NULL;

	frag = &frag_obj->frag;

	if (splat_levels_init(frag, &levels, levels_obj))
		return NULL;

	all_floats = levels.all_floats;
	all_floats = all_floats && splat_check_all_floats(freq, phase, ratio);

	if (all_floats)
		splat_square_floats(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				    PyFloat_AS_DOUBLE(phase) + origin,
				    PyFloat_AS_DOUBLE(ratio));
	else if (splat_square_signals(frag, levels.obj, freq, phase, ratio,
				      origin))
		return NULL;

	Py_RETURN_NONE;
}

PyDoc_STRVAR(splat_triangle_doc,
"triangle(fragment, levels, frequency, phase=0.0, origin=0.0, ratio=0.5)\n"
"\n"
"Generate a triangle wave with the given ``ratio`` over the entire "
"``fragment``.  The ratio is between the duration of the high and "
"low states.\n");

static PyObject *splat_triangle(PyObject *self, PyObject *args)
{
	Fragment *frag_obj;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *phase = splat_zero;
	double origin = 0.0;
	PyObject *ratio = splat_init_source_ratio;

	struct splat_fragment *frag;
	struct splat_levels levels;
	int all_floats;

	if (!PyArg_ParseTuple(args, "O!OO|OdO", &splat_FragmentType, &frag_obj,
			      &levels_obj, &freq, &phase, &origin, &ratio))
		return NULL;

	frag = &frag_obj->frag;

	if (splat_levels_init(frag, &levels, levels_obj))
		return NULL;

	all_floats = levels.all_floats;
	all_floats = all_floats && splat_check_all_floats(freq, phase, ratio);

	if (all_floats)
		splat_triangle_floats(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				      PyFloat_AS_DOUBLE(phase) + origin,
				      PyFloat_AS_DOUBLE(ratio));
	else if (splat_triangle_signals(frag, levels.obj, freq, phase, ratio,
					origin))
		return NULL;

	Py_RETURN_NONE;
}

PyDoc_STRVAR(splat_overtones_doc,
"overtones(fragment, levels, frequency, overtones, phase=0.0, origin=0.0)\n"
"\n"
"Generate a sum of overtones as pure sine waves with the given fundamental "
"``frequency`` and ``levels``.\n"
"\n"
"The ``overtones`` are described with a list of 3-tuples containing the "
"ratio between the overtone and the fundamental frequency, the phase and "
"linear levels: ``(ratio, phase, levels)``.  All these values can be signals, "
"and the levels can either be a single value for all channels or individual "
"values for each channel.  The generation is performed over the entire "
"fragment.\n");

static PyObject *splat_overtones(PyObject *self, PyObject *args)
{
	enum {
		OT_RATIO = 0,
		OT_PHASE,
		OT_LEVELS,
	};
	Fragment *frag_obj;
	PyObject *levels_obj;
	PyObject *freq;
	PyObject *overtones_obj;
	PyObject *phase = splat_zero;
	double origin = 0.0;

	struct splat_fragment *frag;
	struct splat_levels levels;
	struct splat_overtone *overtones;
	struct splat_overtone *ot;
	Py_ssize_t n;
	Py_ssize_t pos;
	int all_floats;
	int ot_all_floats;
	double fl_period;
	int stat = 0;

	if (!PyArg_ParseTuple(args, "O!OOO!|Od", &splat_FragmentType, &frag_obj,
			      &levels_obj, &freq, &PyList_Type, &overtones_obj,
			      &phase, &origin))
		return NULL;

	frag = &frag_obj->frag;

	if (splat_levels_init(frag, &levels, levels_obj))
		return NULL;

	n = PyList_GET_SIZE(overtones_obj);
	overtones = PyMem_Malloc(n * sizeof(struct splat_overtone));

	if (overtones == NULL)
		return PyErr_NoMemory();

	if (PyFloat_Check(freq)) {
		fl_period = 1.0 / PyFloat_AS_DOUBLE(freq);
		origin = fmod(origin, fl_period);
	} else {
		fl_period = 0.0;
	}

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

		if (ot_all_floats && PyFloat_Check(ot->ratio)) {
			ot->fl_ratio = PyFloat_AS_DOUBLE(ot->ratio);
#ifdef SPLAT_FAST
			ot->fl_ratioq = sf_set(ot->fl_ratio);
#endif
		} else {
			ot_all_floats = 0;
		}

		ot->phase = PyTuple_GET_ITEM(ot_params, OT_PHASE);

		if (ot_all_floats && PyFloat_Check(ot->phase)) {
			ot->fl_phase = PyFloat_AS_DOUBLE(ot->phase);

			if (fl_period)
				ot->fl_phase = fmod(ot->fl_phase,
						    fl_period / ot->fl_ratio);

#ifdef SPLAT_FAST
			ot->fl_phaseq = sf_set(ot->fl_phase);
#endif
		} else {
			ot_all_floats = 0;
		}

		ot_levels = PyTuple_GET_ITEM(ot_params, OT_LEVELS);

		if (splat_levels_init(frag, &ot->levels, ot_levels))
			goto free_overtones;

		ot_all_floats = ot_all_floats && ot->levels.all_floats;
	}

	all_floats = all_floats && ot_all_floats;

	if (all_floats)
		splat_overtones_float(frag, levels.fl, PyFloat_AS_DOUBLE(freq),
				      PyFloat_AS_DOUBLE(phase) + origin,
				      overtones, n);
	else if (ot_all_floats)
		stat = splat_overtones_mixed(frag, levels.obj, freq, phase,
					     overtones, n, origin);
	else
		stat = splat_overtones_signal(frag, levels.obj, freq, phase,
					      overtones, n, origin);
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
"This filter applies a decreasing linear envelope over the ``fragment`` with "
" ``k`` and ``p`` arguments as follows, for a sound signal ``s`` "
"at index ``i``:\n"
"\n"
".. math::\n"
"\n"
"   s[i] = \\frac{s[i]}{(1 + \\frac{i}{k})^p}\n"
"\n");

static PyObject *splat_dec_envelope(PyObject *self, PyObject *args)
{
	Fragment *frag_obj;
	double k = 1.0;
	double p = 1.0;

	if (!PyArg_ParseTuple(args, "O!|dd", &splat_FragmentType, &frag_obj,
			      &k, &p))
		return NULL;

	if (k == 0.0) {
		PyErr_SetString(PyExc_ValueError, "k must not be 0");
		return NULL;
	}

	splat_filter_dec_envelope(&frag_obj->frag, k, p);

	Py_RETURN_NONE;
}

PyDoc_STRVAR(splat_reverse_doc,
"reverse(fragment)\n"
"\n"
"Reverse the order of all the samples in ``fragment``.\n");

static PyObject *splat_reverse(PyObject *self, PyObject *args)
{
	Fragment *frag_obj;

	if (!PyArg_ParseTuple(args, "O!", &splat_FragmentType, &frag_obj))
		return NULL;

	splat_filter_reverse(&frag_obj->frag);

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
	Fragment *frag_obj;
	PyObject *delays_list;
	double time_factor = 0.2;
	double gain_factor = 6.0;
	unsigned int seed = 0;

	struct splat_fragment *frag;
	struct splat_delay *delays[SPLAT_MAX_CHANNELS];
	Py_ssize_t n_delays;
	size_t delays_size;
	size_t max_delay;
	size_t max_index;
	size_t d;
	unsigned c;

	if (!PyArg_ParseTuple(args, "O!O!|ddI", &splat_FragmentType, &frag_obj,
			      &PyList_Type, &delays_list, &time_factor,
			      &gain_factor, &seed))
		return NULL;

	if (!seed)
		seed = time(0);

	srand(seed);

	frag = &frag_obj->frag;
	n_delays = PyList_GET_SIZE(delays_list);
	delays_size = n_delays * sizeof(struct splat_delay);

	for (c = 0; c < frag->n_channels; ++c) {
		delays[c] = PyMem_Malloc(delays_size);

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

			delays[c][d].time = c_time * frag->rate;

			if (delays[c][d].time > max_delay)
				max_delay = delays[c][d].time;
		}

		gain = PyFloat_AsDouble(PyTuple_GetItem(pair, 1));

		for (c = 0; c < frag->n_channels; ++c) {
			const double c_gain_dB =
				(gain - gain_factor
				 + (rand() * gain_factor * 2.0 / RAND_MAX));

			delays[c][d].gain = dB2lin(c_gain_dB);
		}
	}

	max_index = frag->length + 1;

	if (splat_frag_grow(frag, (frag->length + max_delay)))
		return NULL;

	splat_filter_reverb(&frag_obj->frag, delays, n_delays, max_index);

	for (c = 0; c < frag->n_channels; ++c)
		PyMem_Free(delays[c]);

	Py_RETURN_NONE;
}

static PyObject *splat_poly_value(PyObject *self, PyObject *args)
{
	PyObject *coefs;
	double x;
	PyObject *db_obj;

	int db;

	if (!PyArg_ParseTuple(args, "O!dO!", &PyTuple_Type, &coefs, &x,
			      &PyBool_Type, &db_obj))
		return NULL;

	db = (db_obj == Py_True) ? 1 : 0;

	return PyFloat_FromDouble(splat_spline_tuple_value(coefs, x, db));
}

static PyObject *splat_spline_value(PyObject *self, PyObject *args)
{
	PyObject *spline;
	double x;
	PyObject *db_obj;

	PyObject *poly;
	int db;

	if (!PyArg_ParseTuple(args, "O!dO!", &PyList_Type, &spline, &x,
			      &PyBool_Type, &db_obj))
		return NULL;

	poly = splat_spline_find_poly(spline, x, NULL);

	if (poly == NULL)
		Py_RETURN_NONE;

	db = (db_obj == Py_True) ? 1 : 0;

	return PyFloat_FromDouble(splat_spline_tuple_value(poly, x, db));
}

static PyMethodDef splat_methods[] = {
	{ "lin2dB", splat_lin2dB, METH_VARARGS,
	  splat_lin2dB_doc },
	{ "dB2lin", splat_dB2lin, METH_VARARGS,
	  splat_dB2lin_doc },
	{ "gen_ref", splat_gen_ref, METH_VARARGS,
	  splat_gen_ref_doc },
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
	{ "poly_value", splat_poly_value, METH_VARARGS, NULL },
	{ "spline_value", splat_spline_value, METH_VARARGS, NULL },
	{ NULL, NULL, 0, NULL }
};

static void splat_init_sample_types(PyObject *m, const char *name,
				    PyObject *obj)
{
	size_t i;

	obj = PyDict_New();

	for (i = 0; i < ARRAY_SIZE(splat_raw_io_table); ++i) {
		const struct splat_raw_io *io = &splat_raw_io_table[i];

		PyDict_SetItem(obj, PyString_FromString(io->sample_type),
			       PyLong_FromLong(io->sample_width));
	}

	PyModule_AddObject(m, name, obj);
}

static void splat_init_page_size(PyObject *m)
{
	Fragment *zero_frag;

	zero_frag = (Fragment *)Fragment_new(&splat_FragmentType, NULL, NULL);

	if (zero_frag == NULL) {
		PyErr_NoMemory();
		return;
	}

	splat_page_size = sysconf(_SC_PAGESIZE);
	if (splat_frag_init(&zero_frag->frag, 1, 8000,
			    splat_page_size / sizeof(sample_t), NULL))
		return;

	splat_zero_page = zero_frag->frag.channels[0].data;
	PyModule_AddObject(m, "_zero_page_frag", (PyObject *)zero_frag);
}

PyMODINIT_FUNC init_splat(void)
{
	struct splat_type {
		PyTypeObject *type;
		const char *name;
	};
	static const struct splat_type splat_types[] = {
		{ &splat_SplineType, "Spline" },
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
	splat_zero = PyFloat_FromDouble(0.0);
	PyModule_AddObject(m, "_zero", splat_zero);
	splat_one = PyFloat_FromDouble(1.0);
	PyModule_AddObject(m, "_one", splat_one);
	splat_init_sample_types(m, "sample_types", splat_sample_types);
	splat_init_page_size(m);

	PyModule_AddStringConstant(m, "SAMPLE_TYPE", SPLAT_NATIVE_SAMPLE_TYPE);
	PyModule_AddIntConstant(m, "SAMPLE_WIDTH", SPLAT_NATIVE_SAMPLE_WIDTH);

#ifdef SPLAT_FAST
	splat_sine_step = sf_set((float)splat_sine_table_len / M_PI);
#endif
}

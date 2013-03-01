#include <Python.h>
#include <math.h>

#define BASE_TYPE_FLAGS (Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE)
#define MAX_CHANNELS 16

#define lin2dB(level) (20 * log10(level))
#define dB2lin(dB) (pow10((dB) / 20))

/* ----------------------------------------------------------------------------
 * Fragment implementation
 */

struct Fragment_object {
	PyObject_HEAD;
	int init;
	unsigned n_channels;
	unsigned rate;
	Py_ssize_t length;
	float *data[MAX_CHANNELS];
};
typedef struct Fragment_object Fragment;

static PyTypeObject geomusic_FragmentType;

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
	unsigned n_channels;
	unsigned rate;
	float duration;

	unsigned i;
	size_t length;
	size_t data_size;

	if (!PyArg_ParseTuple(args, "IIf", &n_channels, &rate, &duration))
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
	data_size = length * sizeof(float);

	for (i = 0; i < n_channels; ++i) {
		self->data[i] = PyMem_Malloc(data_size);

		if (!self->data[i]) {
			PyErr_NoMemory();
			return -1;
		}

		memset(self->data[i], 0, data_size);
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

	if (!self)
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

	if (!sample)
		return PyErr_NoMemory();

	for (c = 0; c < self->n_channels; ++c) {
		const double s = (double)self->data[c][i];
		PyTuple_SetItem(sample, c, PyFloat_FromDouble(s));
	}

	return sample;
}

static int Fragment_sq_ass_item(Fragment *self, Py_ssize_t i, PyObject *v)
{
	Py_ssize_t c;

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

		self->data[c][i] = (float)PyFloat_AS_DOUBLE(s);
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

static PyObject *Fragment_get_sample_rate(Fragment *self, void *_)
{
	return Py_BuildValue("I", self->rate);
}

static PyObject *Fragment_get_duration(Fragment *self, void *_)
{
	return Py_BuildValue("f", (float)self->length / self->rate);
}

static PyObject *Fragment_get_channels(Fragment *self, void *_)
{
	return Py_BuildValue("I", self->n_channels);
}

static PyGetSetDef Fragment_getsetters[] = {
	{ "sample_rate", (getter)Fragment_get_sample_rate, NULL, NULL },
	{ "duration", (getter)Fragment_get_duration, NULL, NULL },
	{ "channels", (getter)Fragment_get_channels, NULL, NULL },
	{ NULL }
};

/* Fragment methods */

static int do_resize(Fragment *self, size_t length)
{
	size_t start;
	size_t data_size;
	size_t i;

	if (length <= self->length)
		return 0;

	start = self->length * sizeof(float);
	data_size = length * sizeof(float);

	for (i = 0; i < self->n_channels; ++i) {
		self->data[i] = PyMem_Realloc(self->data[i], data_size);

		if (!self->data[i]) {
			PyErr_NoMemory();
			return -1;
		}

		memset(&self->data[i][self->length], 0, (data_size - start));
	}

	self->length = length;

	return 0;
}

static PyObject *Fragment_resize(Fragment *self, PyObject *args)
{
	Py_ssize_t length;

	if (!PyArg_ParseTuple(args, "n", &length))
		return NULL;

	if (do_resize(self, length) < 0)
		return NULL;

	Py_RETURN_NONE;
}

static PyObject *Fragment_mix(Fragment *self, PyObject *args)
{
	Fragment *frag;
	float start = 0.0;

	size_t start_sample;
	size_t total_length;
	unsigned c;

	if (!PyArg_ParseTuple(args, "O!|f", &geomusic_FragmentType, &frag,
			      &start))
		return NULL;

	if (frag->n_channels != self->n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return NULL;
	}

	if (frag->rate != self->rate) {
		PyErr_SetString(PyExc_ValueError, "sample rate mismatch");
		return NULL;
	}

	start_sample = start * self->rate;
	total_length = start_sample + frag->length;

	if (do_resize(self, total_length) < 0)
		return NULL;

	for (c = 0; c < self->n_channels; ++c) {
		const float *src = frag->data[c];
		float *dst =  &self->data[c][start_sample];
		Py_ssize_t i = frag->length;

		while (i--)
			*dst++ += *src++;
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
	unsigned bytes_per_sample;
	unsigned n_samples;
	unsigned end;
	unsigned ch;

	if (!PyArg_ParseTuple(args, "O!iIII", &PyByteArray_Type, &bytes_obj,
			      &start, &sample_width, &sample_rate,
			      &n_channels))
		return NULL;

	if (sample_width != 2) {
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

	if (do_resize(self, end) < 0)
		return NULL;

	for (ch = 0; ch < self->n_channels; ++ch) {
		const void *in = bytes + (sample_width * ch);
		float *out = &self->data[ch][start];
		unsigned s;

		for (s = start; s < end; ++s) {
			*out++ = *(int16_t *)in / 32678.0;
			in += bytes_per_sample;
		}
	}

	Py_RETURN_NONE;
}

static PyObject *Fragment_as_bytes(Fragment *self, PyObject *args)
{
	unsigned sample_width;

	PyObject *bytes_obj;
	Py_ssize_t bytes_size;
	const float *it[MAX_CHANNELS];
	unsigned i;
	unsigned c;
	char *out;

	if (!PyArg_ParseTuple(args, "I", &sample_width))
		return NULL;

	if (sample_width != 2) {
		PyErr_SetString(PyExc_ValueError, "unsupported sample width");
		return NULL;
	}

	bytes_obj = PyByteArray_FromStringAndSize("", 0);

	if (!bytes_obj)
		return PyErr_NoMemory();

	bytes_size = self->length * self->n_channels * sample_width;

	if (PyByteArray_Resize(bytes_obj, bytes_size))
		return PyErr_NoMemory();

	for (c = 0; c < self->n_channels; ++c)
		it[c] = self->data[c];

	out = PyByteArray_AS_STRING(bytes_obj);

	for (i = 0; i < self->length; ++i) {
		for (c = 0; c < self->n_channels; ++c) {
			const float z = *(it[c]++);
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

	return bytes_obj;
}

static PyMethodDef Fragment_methods[] = {
	{ "_resize", (PyCFunction)Fragment_resize, METH_VARARGS,
	  "Resize the internal buffer" },
	{ "mix", (PyCFunction)Fragment_mix, METH_VARARGS,
	  "Mix two fragments together" },
	{ "import_bytes", (PyCFunction)Fragment_import_bytes, METH_VARARGS,
	  "Import data as raw bytes" },
	{ "as_bytes", (PyCFunction)Fragment_as_bytes, METH_VARARGS,
	  "Make a byte buffer with the data" },
	{ NULL }
};

static PyTypeObject geomusic_FragmentType = {
	PyObject_HEAD_INIT(NULL)
	0,                                 /* ob_size */
	"_geomusic.Fragment",              /* tp_name */
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
	"Geomusic Fragment",               /* tp_doc */
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
 * _geomusic methods
 */

static PyObject *geomusic_lin2dB(PyObject *self, PyObject *args)
{
	float level;

	if (!PyArg_ParseTuple(args, "f", &level))
		return NULL;

	return PyFloat_FromDouble(lin2dB(level));
}

static PyObject *geomusic_dB2lin(PyObject *self, PyObject *args)
{
	float dB;

	if (!PyArg_ParseTuple(args, "f", &dB))
		return NULL;

	return PyFloat_FromDouble(dB2lin(dB));
}

static PyObject *geomusic_sine(PyObject *self, PyObject *args)
{
	Fragment *frag;
	float freq;
	PyObject *levels_tuple;

	Py_ssize_t n_channels;
	float levels[MAX_CHANNELS];
	Py_ssize_t c, i;
	float k;

	if (!PyArg_ParseTuple(args, "O!fO!", &geomusic_FragmentType, &frag,
			      &freq, &PyTuple_Type, &levels_tuple))
		return NULL;

	n_channels = PyTuple_Size(levels_tuple);

	if (n_channels > MAX_CHANNELS) {
		PyErr_SetString(PyExc_ValueError, "too many channels");
		return NULL;
	}

	if (n_channels != frag->n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return NULL;
	}

	for (c = 0; c < n_channels; ++c) {
		PyObject *level = PyTuple_GetItem(levels_tuple, c);
		levels[c] = (float)dB2lin(PyFloat_AsDouble(level));
	}

	k = 2 * M_PI * freq / frag->rate;

	for (i = 0; i < frag->length; ++i) {
		const float s = sin(k * i);

		for (c = 0; c < n_channels; ++c)
			frag->data[c][i] = s * levels[c];
	}

	Py_RETURN_NONE;
}

static PyObject *geomusic_overtones(PyObject *self, PyObject *args)
{
	struct overtone {
		float freq;
		float levels[MAX_CHANNELS];
	};

	Fragment *frag;
	float freq;
	PyObject *levels_obj;
	PyObject *overtones_obj;

	struct overtone *overtones;
	struct overtone *ot;
	const struct overtone *ot_end;
	Py_ssize_t n_overtones;
	PyObject *ot_freq;
	PyObject *ot_levels;
	float levels[MAX_CHANNELS];
	Py_ssize_t pos;
	size_t i;
	size_t c;
	float k;
	int stat = 0;

	if (!PyArg_ParseTuple(args, "O!fO!O!", &geomusic_FragmentType, &frag,
			      &freq, &PyTuple_Type, &levels_obj,
			      &PyDict_Type, &overtones_obj))
		return NULL;

	if (PyTuple_GET_SIZE(levels_obj) != frag->n_channels) {
		PyErr_SetString(PyExc_ValueError, "channels number mismatch");
		return NULL;
	}

	for (i = 0; i < frag->n_channels; ++i) {
		PyObject *l = PyTuple_GET_ITEM(levels_obj, i);
		levels[i] = (float)dB2lin(PyFloat_AsDouble(l));
	}

	n_overtones = PyDict_Size(overtones_obj);
	overtones = PyMem_Malloc(n_overtones * sizeof(struct overtone));

	if (overtones == NULL)
		return PyErr_NoMemory();

	ot = overtones;
	ot_end = &overtones[n_overtones];
	pos = 0;

	while (PyDict_Next(overtones_obj, &pos, &ot_freq, &ot_levels)) {
		if (!PyFloat_Check(ot_freq)) {
			PyErr_SetString(PyExc_TypeError,
					"overtone key must be a float");
			stat = -1;
			goto free_overtones;
		}

		if (!PyTuple_Check(ot_levels)) {
			PyErr_SetString(PyExc_TypeError,
					"overtone levels must be a tuple");
			stat = -1;
			goto free_overtones;
		}

		if (PyTuple_GET_SIZE(ot_levels) != frag->n_channels) {
			PyErr_SetString(PyExc_ValueError,
					"channels number mismatch");
			stat = -1;
			goto free_overtones;
		}

		ot->freq = (float)PyFloat_AS_DOUBLE(ot_freq);

		for (c = 0; c < frag->n_channels; ++c) {
			PyObject *l = PyTuple_GET_ITEM(ot_levels, c);

			if (!PyFloat_Check(l)) {
				PyErr_SetString(
					PyExc_TypeError,
					"overtone level must be a float");
				stat = -1;
				goto free_overtones;
			}

			ot->levels[c] = ((float)dB2lin(PyFloat_AS_DOUBLE(l))
					 * levels[c]);
		}

		++ot;
	}

	k = 2 * M_PI * freq / frag->rate;

	/* Silence harmonics above rate / 2 to avoid spectrum overlap */
	for (ot = overtones; ot != ot_end; ++ot)
		if ((ot->freq * freq) >= (frag->rate / 2))
			for (c = 0; c < frag->n_channels; ++c)
				ot->levels[c] = 0.0f;

	for (i = 0; i < frag->length; ++i) {
		const float m = k * i;

		for (ot = overtones; ot != ot_end; ++ot) {
			const float s = sin(m * ot->freq);

			for (c = 0; c < frag->n_channels; ++c)
				frag->data[c][i] += s * ot->levels[c];
		}
	}

free_overtones:
	PyMem_Free(overtones);

	if (stat < 0)
		return NULL;

	Py_RETURN_NONE;
}

static PyObject *geomusic_normalize(PyObject *self, PyObject *args)
{
	Fragment *frag;
	float level;

	float average[MAX_CHANNELS];
	unsigned c;
	float peak;
	float gain;

	if (!PyArg_ParseTuple(args, "O!f", &geomusic_FragmentType, &frag,
			      &level))
		return NULL;

	if (frag->n_channels > MAX_CHANNELS) {
		PyErr_SetString(PyExc_ValueError, "too many channels");
		return NULL;
	}

	level = dB2lin(level);
	peak = 0.0;

	for (c = 0; c < frag->n_channels; ++c) {
		float * const chan_data = frag->data[c];
		const float * const end = &chan_data[frag->length];
		const float *it;
		float avg = 0.0;
		float neg = 1.0;
		float pos = -1.0;
		float chan_peak;

		for (it = frag->data[c]; it != end; ++it) {
			avg += *it / frag->length;

			if (*it > pos)
				pos = *it;

			if (*it < neg)
				neg = *it;
		}

		average[c] = avg;
		neg -= avg;
		neg = fabsf(neg);
		pos -= avg;
		pos = fabsf(pos);
		chan_peak = (neg > pos) ? neg : pos;

		if (chan_peak > peak)
			peak = chan_peak;
	}

	gain = level / peak;

	for (c = 0; c < frag->n_channels; ++c) {
		const float chan_avg = average[c];
		float * const chan_data = frag->data[c];
		const float * const end = &chan_data[frag->length];
		float *it;

		for (it = chan_data; it != end; ++it) {
			*it -= chan_avg;
			*it *= gain;
		}
	}

	Py_RETURN_NONE;
}

static PyMethodDef geomusic_methods[] = {
	{ "lin2dB", geomusic_lin2dB, METH_VARARGS,
	  "Convert linear value to dB" },
	{ "dB2lin", geomusic_dB2lin, METH_VARARGS,
	  "Convert dB value to linear" },
	{ "sine", geomusic_sine, METH_VARARGS, "Make a sine wave" },
	{ "overtones", geomusic_overtones, METH_VARARGS,
	  "Add a series of overtones" },
	{ "normalize", geomusic_normalize, METH_VARARGS, "Normalize" },
	{ NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC init_geomusic(void)
{
	struct geomusic_type {
		PyTypeObject *type;
		const char *name;
	};
	static const struct geomusic_type geomusic_types[] = {
		{ &geomusic_FragmentType, "Fragment" },
		{ NULL, NULL }
	};
	const struct geomusic_type *it;
	PyObject *m;

	for (it = geomusic_types; it->type != NULL; ++it) {
		if (it->type->tp_new == NULL)
			it->type->tp_new = PyType_GenericNew;

		if (PyType_Ready(it->type))
			return;
	}

	m = Py_InitModule("_geomusic", geomusic_methods);

	for (it = geomusic_types; it->type != NULL; ++it) {
		Py_INCREF((PyObject *)it->type);
		PyModule_AddObject(m, it->name, (PyObject *)it->type);
	}
}

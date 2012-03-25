#include <Python.h>
#include <math.h>

#define MAX_CHANNELS 16

static PyObject *geomusic_sine(PyObject *self, PyObject *args)
{
	PyObject *data;
	float freq;
	unsigned long rate;
	Py_ssize_t length;
	PyObject *levels_tuple;

	Py_ssize_t n_channels;
	float levels[MAX_CHANNELS];
	PyObject *channels[MAX_CHANNELS];
	Py_ssize_t c, i;
	float k;

	if (!PyArg_ParseTuple(args, "OfknO!", &data, &freq, &rate, &length,
			      &PyTuple_Type, &levels_tuple))
		return NULL;

	n_channels = PyTuple_Size(levels_tuple);

	if (n_channels > MAX_CHANNELS) {
		fprintf(stderr, "Too many channels: %zi\n", n_channels);
		return NULL;
	}

	/* ToDo: levels in dB */
	for (c = 0; c < n_channels; ++c) {
		PyObject *level = PyTuple_GetItem(levels_tuple, c);
		levels[c] = (float)PyFloat_AsDouble(level);
		channels[c] = PyList_GetItem(data, c);
	}

	k = 2 * M_PI * freq / rate;

	for (i = 0; i < length; ++i) {
		const double s = sin(k * i);

		for (c = 0; c < n_channels; ++c) {
			const double z = s * levels[c];
			PyList_SET_ITEM(channels[c], i, PyFloat_FromDouble(z));
		}
	}

	Py_RETURN_NONE;
}

static PyMethodDef geomusic_methods[] = {
	{ "sine", geomusic_sine, METH_VARARGS, "Make a sine wave" },
	{ NULL, NULL, 0, NULL }
};

PyMODINIT_FUNC init_geomusic(void)
{
	PyObject *m;

	m = Py_InitModule("_geomusic", geomusic_methods);
}

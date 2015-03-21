/*
    Splat - spline.c

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

double splat_spline_tuple_value(PyObject *poly, double x, int db)
{
	Py_ssize_t i;
	double value = 0.0;
	double x_pow = 1.0;

	for (i = 0; i < PyTuple_GET_SIZE(poly); ++i) {
		PyObject *py_k = PyTuple_GET_ITEM(poly, i);
		const double k = PyFloat_AS_DOUBLE(py_k);

		value += k * x_pow;
		x_pow *= x;
	}

	if (db)
		value = dB2lin(value);

	return value;
}

PyObject *splat_spline_find_poly(PyObject *spline, double x, double *end)
{
	Py_ssize_t i;

	if (!PyList_CheckExact(spline)) {
		PyErr_SetString(PyExc_TypeError, "spline must be a list");
		return NULL;
	}

	for (i = 0; i < PyList_GET_SIZE(spline); ++i) {
		PyObject *poly_params = PyList_GET_ITEM(spline, i);
		PyObject *param;
		double poly_start;
		double poly_end;

		if (!PyTuple_CheckExact(poly_params) ||
		    (PyTuple_GET_SIZE(poly_params) != 3)) {
			PyErr_SetString(PyExc_TypeError,
					"spline list item must be a 3-tuple");
			return NULL;
		}

		param = PyTuple_GET_ITEM(poly_params, 1);

		if (!PyFloat_CheckExact(param)) {
			PyErr_SetString(PyExc_TypeError,
				"spline list item start time must be a float");
			return NULL;
		}

		poly_end = PyFloat_AS_DOUBLE(param);

		if (x > poly_end)
			continue;

		param = PyTuple_GET_ITEM(poly_params, 0);

		if (!PyFloat_CheckExact(param)) {
			PyErr_SetString(PyExc_TypeError,
				"spline list item end time must be a float");
			return NULL;
		}

		poly_start = PyFloat_AS_DOUBLE(param);

		if (x < poly_start)
			continue;

		param = PyTuple_GET_ITEM(poly_params, 2);

		if (!PyTuple_CheckExact(param)) {
			PyErr_SetString(PyExc_TypeError,
				"spline list item coefs must be a tuple");
			return NULL;
		}

		if (end != NULL)
			*end = poly_end;

		return param;
	}

	return NULL;
}

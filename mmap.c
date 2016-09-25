/*
    Splat - mmap.c

    Copyright (C) 2016, 2017
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

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include "_splat.h"

static int splat_mmap(struct splat_mmap *m, size_t sz)
{
	m->ptr = mmap(NULL, sz, PROT_READ | PROT_WRITE, MAP_SHARED, m->fd, 0);

	if (m->ptr == MAP_FAILED)
		return -1;

	m->size = sz;

	return 0;
}

static void splat_unmap(struct splat_mmap *m)
{
	munmap(m->ptr, m->size);
	m->ptr = NULL;
}

static size_t splat_mmap_n_pages(size_t size)
{
	size_t n_pages;

	n_pages = size / splat_page_size;

	if (!size || (size % splat_page_size))
		n_pages++;

	return n_pages;
}

int splat_mmap_init(struct splat_mmap *m, const char *path)
{
	static const mode_t mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;
	struct stat s;
	int flags = O_RDWR;

	m->ptr = NULL;
	m->size = 0;
	m->persist = 0;

	if ((stat(path, &s) < 0) && (errno == ENOENT))
		flags |= O_CREAT;

	m->fd = open(path, flags, mode);

	if (m->fd < 0)
		return -1;

	m->path = strdup(path);

	if (m->path == NULL) {
		close(m->fd);
		return -1;
	}

	if (!(flags & O_CREAT))
		return splat_mmap(m, lseek(m->fd, 0, SEEK_END));

	return 0;
}

void splat_mmap_free(struct splat_mmap *m)
{
	if (m->ptr)
		munmap(m->ptr, m->size);

	close(m->fd);

	if (!m->persist)
		unlink(m->path);

	free(m->path);
}

static int splat_mmap_remap_new(struct splat_mmap *m, size_t sz, size_t pages)
{
	if (lseek(m->fd, 0, SEEK_END) != 0)
		return -1;

	while (pages--) {
		if (write(m->fd, splat_zero_page, splat_page_size)
		    != splat_page_size)
			return -1;
	}

	return splat_mmap(m, sz);
}

static int splat_mmap_remap_grow(struct splat_mmap *m, size_t sz, size_t pages)
{
	const size_t current_pages = m->size / splat_page_size;
	size_t extra_pages = pages - current_pages;

	if (!extra_pages)
		return 0;

	splat_unmap(m);

	if (lseek(m->fd, 0, SEEK_END) != m->size)
		return -1;

	while (extra_pages--) {
		if (write(m->fd, splat_zero_page, splat_page_size)
		    != splat_page_size)
			return -1;
	}

	return splat_mmap(m, sz);
}

static int splat_mmap_remap_shrink(struct splat_mmap *m, size_t sz,
				   size_t length)
{
	memset(m->ptr + length, 0, sz - length);

	if (sz == m->size)
		return 0;

	splat_unmap(m);

	if (ftruncate(m->fd, sz) < 0)
		return -1;

	return splat_mmap(m, sz);
}

int splat_mmap_remap(struct splat_mmap *m, size_t length)
{
	const size_t pages = splat_mmap_n_pages(length);
	const size_t sz = pages * splat_page_size;
	const off_t file_size = lseek(m->fd, 0, SEEK_END);
	int stat;

	if (!file_size) {
		stat = splat_mmap_remap_new(m, sz, pages);
	} else if (sz > m->size) {
		stat = splat_mmap_remap_grow(m, sz, pages);
	} else {
		stat = splat_mmap_remap_shrink(m, sz, length);
	}

	return stat;
}

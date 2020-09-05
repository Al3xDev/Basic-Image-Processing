/*
Copyright 2020 Stoica Alexandru-Gabriel

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

*/


#ifndef BITMAP_H
#define BITMAP_H

#include <stdint.h>

typedef union pixel
{
    uint8_t data[3];
    struct
    {
        uint8_t b,g,r;
    };
} pixel;

typedef struct bitmap
{
    int32_t width,height;
    pixel *data;

} bitmap;

bitmap* bitmap_load( char*);
void bitmap_save( bitmap*, char*);

void bitmap_get_pixel( bitmap*, int, int, pixel*);
void bitmap_set_pixel( bitmap*, int, int, pixel*);

#endif

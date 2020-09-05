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

#include "bitmap.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

typedef struct// it must be 14 bytes
{               // real size is 16
    uint16_t sig_code;
    uint32_t file_size;
    uint32_t res;
    uint32_t data_offset;

} bitmap_header;

typedef struct
{
    uint32_t DIB_header_size;
    int32_t image_width;
    int32_t image_height;
    uint16_t planes;
    uint16_t bits_per_pixel;
    uint32_t compresion;
    uint32_t image_size;
    int32_t x_pixels_per_meter;
    int32_t y_pixels_per_meter;
    uint32_t colors_used;
    uint32_t colors_important;

} image_header;


bitmap* bitmap_load( char *name)
{
    FILE *file = fopen( name, "rb");

    bitmap_header bh;
    image_header ih;
    bitmap *bmp;

    fread( (void*)&( bh.sig_code), 1, 2, file);
    fread( (void*)&( bh.file_size), 1, 4, file);
    fread( (void*)&( bh.res), 1, 4, file);
    fread( (void*)&( bh.data_offset), 1, 4, file);

    if( bh.sig_code != 0x4d42 && bh.file_size ==0)
    {
        bmp = 0;
        return bmp;
    }

    fread( (void*)&(ih.DIB_header_size), 1, 4,file);
    fread( (void*)&(ih.image_width), 1, 4,file);
    fread( (void*)&(ih.image_height), 1, 4,file);
    fread( (void*)&(ih.planes), 1, 2,file);
    fread( (void*)&(ih.bits_per_pixel), 1, 2,file);
    fread( (void*)&(ih.compresion), 1, 4,file);
    fread( (void*)&(ih.image_size), 1, 4,file);
    fread( (void*)&(ih.x_pixels_per_meter), 1, 4,file);
    fread( (void*)&(ih.y_pixels_per_meter), 1, 4,file);
    fread( (void*)&(ih.colors_used), 1, 4,file);
    fread( (void*)&(ih.colors_important), 1, 4,file);

    if( ih.bits_per_pixel != 24 && ih.compresion != 0 && ih.image_height == 0 && ih.image_width == 0)
    {
        bmp = 0;
        return bmp;
    }

    bmp = ( bitmap*)malloc( sizeof( bitmap));
    bmp->width = ih.image_width;
    bmp->height = ih.image_height;
    bmp->data = ( pixel*) malloc( sizeof( pixel) * bmp->width * bmp->height);
    fseek( file, bh.data_offset-4, SEEK_SET);

    int padding =  (sizeof( pixel) * bmp->width * bmp->height) % 4;
    if( padding != 0) padding = 4 - padding;
    int padding_buffer;


    for( int row = 0; row < bmp->height; ++row)
    {
        fread( (void*)( bmp->data + bmp->width * row), 1, sizeof( pixel)* bmp->width, file);
        if( padding !=0 ) fread( (void*)&padding_buffer, 1, padding, file);
    }

    return bmp;
}


void bitmap_save( bitmap *bmp,char *name)
{
    int padding_buffer = 0;
    int padding =  (sizeof( pixel) * bmp->width * bmp->height)%4;
    if( padding != 0) padding = 4 - padding;


    bitmap_header bh;
    image_header ih;

    bh.sig_code = 0x4d42;
    bh.file_size = 14 + sizeof( image_header) + sizeof( pixel)* bmp->width * bmp->height + bmp->height*padding;
    bh.data_offset = 0 + 14 + sizeof( image_header);

    ih.DIB_header_size = 40;
    ih.image_width = bmp->width;
    ih.image_height = bmp->height;
    ih.planes = 1;
    ih.bits_per_pixel = 24;
    ih.compresion = 0;
    ih.image_size = bmp->width * bmp->height * sizeof(pixel);
    ih.x_pixels_per_meter = 0;
    ih.y_pixels_per_meter = 0;
    ih.colors_used = 256*256*256;
    ih.colors_important = 256*256*256;


    FILE *file = fopen( name, "wr");

    fwrite( (void*)&bh.sig_code, 1, 2, file);
    fwrite( (void*)&bh.file_size, 1, 4, file);
    fwrite( (void*)&bh.res, 1, 4, file);
    fwrite( (void*)&bh.data_offset, 1, 4, file);

    fwrite( (void*)&(ih.DIB_header_size), 1, 4,file);
    fwrite( (void*)&(ih.image_width), 1, 4,file);
    fwrite( (void*)&(ih.image_height), 1, 4,file);
    fwrite( (void*)&(ih.planes), 1, 2,file);
    fwrite( (void*)&(ih.bits_per_pixel), 1, 2,file);
    fwrite( (void*)&(ih.compresion), 1, 4,file);
    fwrite( (void*)&(ih.x_pixels_per_meter), 1, 4,file);
    fwrite( (void*)&(ih.y_pixels_per_meter), 1, 4,file);
    fwrite( (void*)&(ih.colors_used), 1, 4,file);
    fwrite( (void*)&(ih.colors_important), 1, 4,file);


    for( int row =0; row < bmp->height; ++row)
    {

        fwrite( (void*)(bmp->data + bmp->width*row), 1, sizeof(pixel) * bmp->width, file);
        if( padding != 0) fwrite( (void*) &padding_buffer, 1, padding, file);
    }

}


void bitmap_get_pixel( bitmap *bmp, int x, int y, pixel *p)
{
    p->r = bmp->data[ x + y*bmp->width].r;
    p->g = bmp->data[ x + y*bmp->width].g;
    p->b = bmp->data[ x + y*bmp->width].b;
}

void bitmap_set_pixel( bitmap *bmp, int x, int y, pixel *p)
{
    bmp->data[ x + y*bmp->width].r = p->r;
    bmp->data[ x + y*bmp->width].g = p->g;
    bmp->data[ x + y*bmp->width].b = p->b;
}

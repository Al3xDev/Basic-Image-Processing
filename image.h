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

#ifndef IMAGE_H
#define IMAGE_H
#include <cstdint>

class img_grey
{
    public:
        img_grey();
        void create( int32_t width, int32_t height);
        img_grey* arthmetic_mean_filer( int32_t N);
        img_grey* geometric_mean_filer( int32_t N);
        float* histrogram();
        void histrogram_equ();
        void histrogram_equ_by_part(int);
        void contrast( float);
        void auto_scale();/** rezultat slabut **/
        void auto_scale_ragion( int, int, int, int);
        void add( img_grey*,float);
        void brig( int32_t);
        void destroy();

        img_grey* app_kernel( float*, int);

    int32_t width,height;
    uint8_t *data;
};

struct pixel_t
{
    uint8_t r,g,b;

};

#define R 0
#define G 1
#define B 2

class img_rgb
{
    public:
        img_rgb();
        img_rgb* create( int32_t width, int32_t height);
        img_rgb* geometric_mean_filter( int32_t N);
        float* histrogram();

        void destroy();

    int32_t width,height;
    pixel_t *data;
};




#endif // IMAGE_H

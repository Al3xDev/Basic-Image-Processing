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

#include "image.h"
#include <cmath>
#include <ctime>
#include <iostream>
/*IMAGE_GREY */

img_grey::img_grey() { width=0; height=0; data=0;}

img_grey *img_grey::app_kernel( float *kernel, int size)
{
    if(size%2==0)
    {
        ///error
        return 0;
    }
    img_grey *out;
    out = new img_grey;
    out->height = height;
    out->width = width;
    out->data = new uint8_t[width*height];


    uint8_t *new_data = out->data;

    int size_2 = size/2;
    for( int i=size_2; i<width-size_2; ++i)
    {
        for( int j=size_2; j<height-size_2; ++j)
        {
            new_data[ i + j*width] = 0;
            for( int q=-size/2; q<=size/2; ++q)
            {
                for( int w=-size/2; w<=size/2; ++w)
                {
                    new_data[ i + j*width] += data[ (i+w)+(j+q)*width]*kernel[ size/2+w + (q+size/2)*size];
                }
            }
        }
    }

    return out;
}

void img_grey::create( int32_t _width, int32_t _height)
{
    if(_width<=0 || _height<=0)
    {
        /* error msg*/
        return ;
    }

    width = _width;
    height = _height;

    data = new uint8_t[ width*height];
}

img_grey* img_grey::arthmetic_mean_filer( int32_t N)
{
    if( width<=0 || height <=0 || data==0 || ( width < 2*N+1) || ( height < 2*N+1))
    {
        /* error msg*/
        return 0;
    }

    float NN = (2*N+1); NN = NN*NN;

    img_grey *img = new img_grey;
    img->create( width, height);

    for( int32_t x=0; x<=width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            int32_t sum=0;
            for( int32_t xx=-N; xx<=N; ++xx)
                for( int32_t yy=-N; yy<=N; ++yy)
                    sum += data[ x+xx+(y+yy)*width];
            sum = (float)sum/NN;
            if( N ==0) img->data[ x + y*width] = data[ x + y*width];
            else    img->data[ x + y*width] = sum <=255 ? sum: 255;
        }
    }
    return img;
}


img_grey* img_grey::geometric_mean_filer( int32_t N)
{
    if( width<=0 || height <=0 || data==0 || ( width < 2*N+1) || ( height < 2*N+1))
    {
        /* error msg*/
        return 0;
    }

    float NN = (2*N+1); NN = NN*NN;

    img_grey *img = new img_grey;
    img->create( width, height);

    for( int32_t x=0; x<=width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            float prod=1;
            for( int32_t xx=-N; xx<=N; ++xx)
                for( int32_t yy=-N; yy<=N; ++yy)
                    prod *= pow( data[ x+xx+(y+yy)*width],1.0/NN);
            if( N ==0) img->data[ x + y*width] = data[ x + y*width];
            else    img->data[ x + y*width] = prod <= 255 ? prod: 255;
        }
    }
    return img;
}

void img_grey::add( img_grey *img, float alpha)
{
    if( width!=img->width || height!=img->height || width<=0 || height<=0 || data==0)
    {
        /* error msg */
        std::cerr<<"Error\n\n";
        return;
    }
    float t;
    for( int32_t x=0; x<width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            t=0;
            t = data[ x + y*width] + alpha*img->data[ x + y*width];
            t = t / 2.0;
            if(t > 255)  data[x+y*width]=255;
            else
                data[x+y*width] = t;
        }
    }
}

float* img_grey::histrogram()
{
    float *H =new float[256];
    for( int32_t i=0; i<256;++i) H[i]=0;
    for( int32_t x=0; x<width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            ++H[ data[ x+y*width]];
        }
    }
    for( int32_t i=0; i<256; ++i)
    {
        H[i] = H[i]/(width*height);
    }
    return H;
}


void img_grey::histrogram_equ()
{
    float *H = histrogram();
    float *S_H = new float[256];
    float sum=0;
    for( int i=0; i<256; ++i)
    {
        sum = sum + H[i];
        S_H[ i] = sum;
    }

    double alpha = 256;

    for( int i=0; i<width; ++i)
    {
        for( int j=0; j<height; ++j)
        {
            data[ j*width + i] = alpha * S_H[ data[ j*width  + i]];
        }
    }

    delete[] H;
    delete[] S_H;
}

void img_grey::histrogram_equ_by_part( int n)
{
    float H[256];
    int count;
    int first = true;
    time_t t1,t2;
    for( int k=0; k<256; ++k) H[k]=0;

    int dx,dy;
    dx = width/n;
    dy = height/n;

    /// get histogram n*n
    for( int i=0; i<n; ++i)
    {
        if(first) t1=clock();
        for( int j=0; j<n; ++j)
        {
            count = 0;
            for( int x=0; x<dx; ++x)
            {
                for( int y=0; y<dy; ++y)
                {
                    ++H[ data[(i*dx+x)+(j*dy+y)*width]];
                }
            }

            float S_H[256];
            float sum=0;
            float alpha=256;

            for( int l=0; l<256; ++l)
            {
                sum = sum + H[l];
                S_H[ l] = sum/(n*n);
                H[l] = 0;
            }

            for( int x=0; x<dx; ++x)
            {
                for( int y=0; y<dy; ++y)
                {
                    data[ (i*dx+x) + (j*dy+y)*width] = alpha * S_H[ data[(i*dx+x) + (j*dy+y)*width]];
                }
            }

        }
        if(first)
            {
                t2 = clock();
                float t_for_1 = float(t2-t1)/CLOCKS_PER_SEC;
                std::cerr<<"Time aprox: "<<t_for_1*(width-n)/60<<" min\n";
                first=false;
            }
    }


}

void img_grey::auto_scale()
{
    int32_t min=255;
    int32_t max=0;
    for( int32_t x=0; x<=width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            if( data[ x+y*width] <min) min = data[x+y*width];
            if( data[ x+y*width] >max) max = data[x+y*width];
        }
    }

    for( int32_t x=0; x<width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            data[ x+y*width] = (data[ x+y*width]-min)*255/(max-min);
        }
    }
}

void img_grey::auto_scale_ragion( int _x, int _y, int _width, int _height)
{
    int32_t min=255;
    int32_t max=0;
    for( int32_t x=_x; x<_x+_width; ++x)
    {
        for( int32_t y=_y; y<_y+_height; ++y)
        {
            if( data[ x+y*width] <min) min = data[x+y*width];
            if( data[ x+y*width] >max) max = data[x+y*width];
        }
    }

    if( min!= max)
        for( int32_t x=_x; x<_x+_width; ++x)
        {
            for( int32_t y=_y; y<_y+_height; ++y)
            {
                data[ x+y*width] = (data[ x+y*width]-min)*255/(max-min);
            }
        }


}

void img_grey::brig( int32_t level)
{
    int32_t t;
    for( int32_t x=0; x<width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            t = data[ x+y*width] + level;
            if(t<0) t=0;
            else
                if(t>255) t=255;

            data[ x+y*width] = t;
        }
    }
}

void img_grey::contrast(float alpha)
{
    float avg=0;
    int32_t min = 255;
    for( int32_t x=0; x<=width; ++x)
    {
        int32_t avg_s=0;
        for( int32_t y=0; y<height; ++y)
        {
            avg_s = data[ x+y*width];
            if( data[ x+y*width] <min) min = data[x+y*width];
        }
        avg += (float)avg_s/width/height;
    }

    int32_t t=0;
    for( int32_t x=0; x<=width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            t = (data[ x+y*width]-avg)*alpha+avg;
            if( t>255)
            {
                data[ x+y*width] = 255;
            }
            else
            {
                if( t<0 )
                    data[ x +y*width] = 0;
                else
                    data[ x+y*width] = t;
            }

        }
    }
}

void img_grey::destroy()
{
    if( width<=0 || height<=0)
    {
        /* error msg*/
        return;
    }

    delete[] data;
}


/* IMAGE_RGB */

img_rgb::img_rgb() { width=0; height=0; }

img_rgb* img_rgb::create( int32_t _width, int32_t _height)
{
    if(_width<=0 || _height<=0)
    {
        /* error msg*/
        return 0;
    }

    width = _width;
    height = _height;

    data = new pixel_t[ width*height];
    return this;
}



void img_rgb::destroy()
{
    if( width<=0 || height<=0)
    {
        /* error msg*/
        return;
    }

    delete[] data;
}

float* img_rgb::histrogram()
{
    float *H =new float[ 3*256];

    for( int i=0; i<3*256; ++i) H[i]=0;

    for( int i=0; i<width; ++i)
    {
        for( int j=0; j<height; ++j)
        {
            ++H[ 0*256 + data[ j*width + i].r];
            ++H[ 1*256 + data[ j*width + i].g];
            ++H[ 2*256 + data[ j*width + i].b];
        }
    }

    for( int i=0; i<3*256; ++i)
    {
        H[i] = H[i]/(width*height);
    }
    return H;
}

img_rgb* img_rgb::geometric_mean_filter( int32_t N)
{
    if( width<=0 || height <=0 || data==0 || ( width < 2*N+1) || ( height < 2*N+1))
    {
        /* error msg*/
        return 0;
    }

    float NN = (2*N+1); NN = NN*NN;

    img_rgb *img = new img_rgb;
    img->create( width, height);

    for( int32_t x=0; x<=width; ++x)
    {
        for( int32_t y=0; y<height; ++y)
        {
            float prodr=1;
            float prodg=1;
            float prodb=1;

            for( int32_t xx=-N; xx<=N; ++xx)
                for( int32_t yy=-N; yy<=N; ++yy)
                    {
                        prodr *= pow( data[ x+xx+(y+yy)*width].r,1.0/NN);
                        prodg *= pow( data[ x+xx+(y+yy)*width].g,1.0/NN);
                        prodb *= pow( data[ x+xx+(y+yy)*width].b,1.0/NN);
                    }
            if( N ==0) {
                img->data[ x + y*width].r = data[ x + y*width].r;
                img->data[ x + y*width].g = data[ x + y*width].g;
                img->data[ x + y*width].b = data[ x + y*width].b;
            }
            else
            {
                img->data[ x + y*width].r = prodr <= 255 ? prodr: 255;
                img->data[ x + y*width].g = prodg <= 255 ? prodg: 255;
                img->data[ x + y*width].b = prodb <= 255 ? prodb: 255;
            }
        }
    }
    return img;
}

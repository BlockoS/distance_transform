/*
 Distance_transform computation using the fast sweeping method 
 =========================================================================
 v0.1.1
 Licensed under the MIT License
 (c) 2020-2025 Vincent Cruz
 
 This header file provides a function to compute the distance transform of an image.
 The implementation is based upon the following papers:
 
  1. "A fast sweeping method for Eikonal equations"
     by H. Zhao, Mathematics of computation, 74 (2005),pp. 603–627
     https://www.math.uci.edu/~zhao/homepage/research_files/FSM.pdf
 
  2. "Finding the Skeleton of 2D Shape and Contours: Implementation of Hamilton-Jacobi Skeleton"
     by Yuchen He, Sung Ha Kang, Luis Álvarez (2020)
     http://www.ipol.im/pub/art/2021/296/
  
 * Building:
 Before including this file, add the following line in the file where you want to have the 
 implementation.
	#define DISTANCE_TRANSFORM_IMPLEMENTATION
 
 * Usage:
 
 	void distance_transform(const uint8_t *in, float *out, int width, int height);

 "in" is the pointer an 8bpp greyscale image.
 "output" is the pointer to the a floating point array of at least width*height elements. 
 It will contain the distance map.
 "width" and "height" are its dimension.
 All pixels with a value of 0 are considered to be "inside", whereas any non-zero pixel is
 considered outside. As a consequence, a pixel is considered to be on the boundary if its
 value is 0 and at least one of its 8-neighbour is not zero.
  
 * Note:
 This piece of code is not meant to be "production ready".

 */

#ifndef DISTANCE_TRANSFORM_INCLUDE_H
#define DISTANCE_TRANSFORM_INCLUDE_H

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

void distance_transform(const uint8_t *in, float *output, int width, int height) __attribute__((nonnull(1,2)));

#ifdef __cplusplus
}
#endif

#endif /* DISTANCE_TRANSFORM_INCLUDE_H */

#ifdef DISTANCE_TRANSFORM_IMPLEMENTATION

static inline __attribute__((always_inline)) float* distance_update(float *ptr, float dx, float dy, int inc) {
	float distance;
	// Eq(2.4) in (1)
	// Eq(4) in (2) 
	if(fabsf(dx-dy) < 1) {
		distance = (dx + dy + sqrtf(2.F-(dx-dy)*(dx-dy))) / 2.F;
	} else {
		distance = fminf(dx, dy) + 1.F;
	}
	*ptr = fminf(*ptr, distance);
	return ptr + inc;
}

void distance_transform(const uint8_t *in, float *out, int width, int height) {
	int i, j;
	uint8_t c;
	const float max_dist = (float)(width*width + height*height);
	float *ptr = out;
	float *line;

	// Initialize distance map by setting distance at boundary point to 0 and any other points to the maximal possible distance (width^2 + height^2).
	// A boundary point is a pixel with a value of 0 and with a least one of its 8-neightbours different from 0.
	c = ~in[0] & (in[1] | in[width] | in[width+1]);
	++in;
	*ptr++ = c ? 0.F : max_dist;

	for(i=1; i<(width-1); i++) {
		c = ~in[0] & (in[-1] | in[1] | in[width-1] | in[width] | in[width+1]);
		++in;
		*ptr++ = c ? 0.F : max_dist;
	}

	c = ~in[0] & (in[-1] | in[width-1] | in[width]);
	++in;
	*ptr++ = c ? 0.F : max_dist;

	for(j=1; j<(height-1); j++) {
		c = ~in[0] & (in[-width] | in[-width+1] | in[1] | in[width] | in[width+1]);
		++in;
		*ptr++ = c ? 0.F : max_dist;
	
		for(i=1; i<(width-1); i++) {
			c = ~in[0] & (in[-1-width] | in[-width] | in[-width+1] | in[-1] | in[1] | in[width-1] | in[width] | in[width+1]);
			++in;
			*ptr++ = c ? 0.F : max_dist;
		}
		c = ~in[0] & (in[-1-width] | in[-width] | in[-1] | in[width-1] | in[width]);
		++in;
		*ptr++ = c ? 0.F : max_dist;
	}
	
	c = ~in[0] & (in[-width] | in[-width+1] | in[1]);
	++in;
	*ptr++ = c ? 0.F : max_dist;
	
	for(i=1; i<(width-1); i++) {
		c = ~in[0] & (in[-1-width] | in[-width] | in[-width+1] | in[-1] | in[1]);
		++in;
		*ptr++ = c ? 0.F : max_dist;
	}
	c = ~in[0] & (in[-1-width] | in[-width] | in[-1]);
	++in;
	*ptr++ = c ? 0.F : max_dist;

	// Update distance in the 4 sweeping directions.
	// 1. y=[0,height[ x=[0,width[
	ptr = out;
	ptr = distance_update(ptr, ptr[1], ptr[width], +1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[width], +1);
	}
	ptr = distance_update(ptr, ptr[-1], ptr[width], +1);

	line = ptr;
	for(j=1; j<(height-1); j++, line+=width) {
		ptr = line;
		ptr = distance_update(ptr, ptr[1], fminf(ptr[-width], ptr[width]), +1);
		for(i=1; i<(width-1); i++) {
			ptr = distance_update( ptr,fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), +1);
		}
		ptr = distance_update(ptr, ptr[-1], fminf(ptr[-width], ptr[width]), +1);
	}

	ptr = distance_update(ptr, ptr[1], ptr[-width], +1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[-width], +1);
	}
	ptr = distance_update(ptr,ptr[-1], ptr[-width], +1);

	// 2. y=[0,height[ x=]width,0]
	ptr = out+width-1;
	ptr = distance_update(ptr, ptr[-1], ptr[width], -1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[width], -1);
	}
	ptr = distance_update(ptr, ptr[1], ptr[width], -1);

	line = out+width;
	for(j=1; j<(height-1); j++, line+=width) {
        ptr = line+width-1;

		ptr = distance_update(ptr, ptr[-1], fminf(ptr[-width], ptr[width]), -1);
		for(i=1; i<(width-1); i++) {
			ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), -1);
		}
		ptr = distance_update(ptr, ptr[1], fminf(ptr[-width], ptr[width]), -1);
	}

	ptr = line + width-1;
	ptr = distance_update(ptr, ptr[-1], ptr[-width], -1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[-width], -1);
	}
	ptr = distance_update(ptr, ptr[1], ptr[-width], -1);

	// 3. y=]height,0] x=]width,0]
	ptr = out+width-1+(height-1)*width;
	ptr = distance_update(ptr, ptr[-1], ptr[-width], -1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[-width], -1);
	}
	ptr = distance_update(ptr, ptr[1], ptr[-width], -1);

	for(j=1; j<(height-1); j++) {
		ptr = distance_update(ptr, ptr[-1], fminf(ptr[-width], ptr[width]), -1);
		for(i=1; i<(width-1); i++) {
			ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), -1);
		}
		ptr = distance_update(ptr, ptr[1], fminf(ptr[-width], ptr[width]), -1);
	}

	ptr = distance_update(ptr, ptr[-1], ptr[width], -1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[width], -1);
	}
	ptr = distance_update(ptr, ptr[1], ptr[width], -1);

	// 4. y=]height,0] x=[0,width[
	line = ptr = out+(height-1)*width;
	ptr = distance_update(ptr, ptr[1], ptr[-width], +1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[-width], +1);
	}
	ptr = distance_update(ptr, ptr[-1], ptr[-width], +1);

	line -= width;
	for(j=1; j<(height-1); j++, line-=width) {
		ptr = line;
		ptr = distance_update(ptr, ptr[1], fminf(ptr[-width], ptr[width]), +1);
		for(i=1; i<(width-1); i++) {
			ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), +1);
		}
		ptr = distance_update(ptr, ptr[-1], fminf(ptr[-width], ptr[width]), +1);
	}

	ptr = line;
	ptr = distance_update(ptr, ptr[1], ptr[width], +1);
	for(i=1; i<(width-1); i++) {
		ptr = distance_update(ptr, fminf(ptr[-1], ptr[1]), ptr[width], +1);
	}
	ptr = distance_update(ptr, ptr[-1], ptr[width], +1);
}

#endif /* DISTANCE_TRANSFORM_IMPLEMENTATION */

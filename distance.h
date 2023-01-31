/*
 Distance_transform computation using the fast sweeping method 
 =========================================================================
 v0.1.0
 Licensed under the MIT License
 (c) 2020-2023 Vincent Cruz
 
 This header file provides a function to compute the distance transform of an image.
 The implementation is based upon the following papers:
 
 "A fast sweeping method for Eikonal equations"
 by H. Zhao, Mathematics of computation, 74 (2005),pp. 603–627
 https://www.math.uci.edu/~zhao/homepage/research_files/FSM.pdf
 
 "Finding the Skeleton of 2D Shape and Contours: Implementation of Hamilton-Jacobi Skeleton"
 by Yuchen He, Sung Ha Kang, Luis Álvarez (2020)
 https://www.ipol.im/pub/pre/296/preprint.pdf
  
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

void distance_transform(const uint8_t *in, float *output, int width, int height);

#ifdef __cplusplus
}
#endif

#endif /* DISTANCE_TRANSFORM_INCLUDE_H */

#ifdef DISTANCE_TRANSFORM_IMPLEMENTATION

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
	*ptr++ = c ? 0.f : max_dist;

	for(i=1; i<(width-1); i++) {
		c = ~in[0] & (in[-1] | in[1] | in[width-1] | in[width] | in[width+1]);
		++in;
		*ptr++ = c ? 0.f : max_dist;
	}

	c = ~in[0] & (in[-1] | in[width-1] | in[width]);
	++in;
	*ptr++ = c ? 0.f : max_dist;

	for(j=1; j<(height-1); j++) {
		c = ~in[0] & (in[-width] | in[-width+1] | in[1] | in[width] | in[width+1]);
		++in;
		*ptr++ = c ? 0.f : max_dist;
	
		for(i=1; i<(width-1); i++) {
			c = ~in[0] & (in[-1-width] | in[-width] | in[-width+1] | in[-1] | in[1] | in[width-1] | in[width] | in[width+1]);
			++in;
			*ptr++ = c ? 0.f : max_dist;
		}
		c = ~in[0] & (in[-1-width] | in[-width] | in[-1] | in[width-1] | in[width]);
		++in;
		*ptr++ = c ? 0.f : max_dist;
	}
	
	c = ~in[0] & (in[-width] | in[-width+1] | in[1]);
	++in;
	*ptr++ = c ? 0.f : max_dist;
	
	for(i=1; i<(width-1); i++) {
		c = ~in[0] & (in[-1-width] | in[-width] | in[-width+1] | in[-1] | in[1]);
		++in;
		*ptr++ = c ? 0.f : max_dist;
	}
	c = ~in[0] & (in[-1-width] | in[-width] | in[-1]);
	++in;
	*ptr++ = c ? 0.f : max_dist;

	// Update distance in the 4 sweeping directions.
#define update_distance(dx,dy,inc) \
do { \
	ptr[0] = fminf(ptr[0], (fabs((dx)-(dy)) >= 1) ? (fminf((dx), (dy)) + 1.f) : (((dx) + (dy) + sqrtf(2.f - ((dx)-(dy))*((dx)-(dy)))) / 2.f)); \
	ptr += (inc); \
} while(0)

	// y=[0,height[ x=[0,width[
	ptr = out;
	update_distance(ptr[1], ptr[width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[width], +1);
	}
	update_distance(ptr[-1], ptr[width], +1);

	line = ptr;
	for(j=1; j<(height-1); j++, line+=width) {
		ptr = line;
		update_distance(ptr[1], fminf(ptr[-width], ptr[width]), +1);
		for(i=1; i<(width-1); i++) {
			update_distance(fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), +1);
		}
		update_distance(ptr[-1], fminf(ptr[-width], ptr[width]), +1);
	}

	update_distance(ptr[1], ptr[-width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[-width], +1);
	}
	update_distance(ptr[-1], ptr[-width], +1);

	// y=[0,height[ x=]width,0]
	ptr = out+width-1;
	update_distance(ptr[-1], ptr[width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[width], -1);
	}
	update_distance(ptr[1], ptr[width], -1);

	line = out+width;
	for(j=1; j<(height-1); j++, line+=width) {
        ptr = line+width-1;

		update_distance(ptr[-1], fminf(ptr[-width], ptr[width]), -1);
		for(i=1; i<(width-1); i++) {
			update_distance(fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), -1);
		}
		update_distance(ptr[1], fminf(ptr[-width], ptr[width]), -1);
	}

	ptr = line + width-1;
	update_distance(ptr[-1], ptr[-width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[-width], -1);
	}
	update_distance(ptr[1], ptr[-width], -1);

	// y=]height,0] x=]width,0]
	ptr = out+width-1+(height-1)*width;
	update_distance(ptr[-1], ptr[-width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[-width], -1);
	}
	update_distance(ptr[1], ptr[-width], -1);

	for(j=1; j<(height-1); j++) {
		update_distance(ptr[-1], fminf(ptr[-width], ptr[width]), -1);
		for(i=1; i<(width-1); i++) {
			update_distance(fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), -1);
		}
		update_distance(ptr[1], fminf(ptr[-width], ptr[width]), -1);
	}

	update_distance(ptr[-1], ptr[width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[width], -1);
	}
	update_distance(ptr[1], ptr[width], -1);

	// y=]height,0] x=[0,width[
	line = ptr = out+(height-1)*width;
	update_distance(ptr[1], ptr[-width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[-width], +1);
	}
	update_distance(ptr[-1], ptr[-width], +1);

	line -= width;
	for(j=1; j<(height-1); j++, line-=width) {
		ptr = line;
		update_distance(ptr[1], fminf(ptr[-width], ptr[width]), +1);
		for(i=1; i<(width-1); i++) {
			update_distance(fminf(ptr[-1], ptr[1]), fminf(ptr[-width], ptr[width]), +1);
		}
		update_distance(ptr[-1], fminf(ptr[-width], ptr[width]), +1);
	}

	ptr = line;
	update_distance(ptr[1], ptr[width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fminf(ptr[-1], ptr[1]), ptr[width], +1);
	}
	update_distance(ptr[-1], ptr[width], +1);

#undef update_distance
}

#endif /* DISTANCE_TRANSFORM_IMPLEMENTATION */

/*
 Distance_transform computation using the fast sweeping method 
 =========================================================================
 v0.0.1
 Licensed under the MIT License
 (c) 2020 Vincent Cruz
 
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

 You can define DISTANCE_TRANSFORM_MALLOC and DISTANCE_TRANSFORM_FREE before the include
 to replace malloc and free.
 
 * Usage:
 
 	float* distance_transform(const uint8_t *in, int width, int height);

 "in" is the pointer an 8bpp greyscale image.
 "width" and "height" are its dimension.
 All pixels with a value of 0 are considered to be "inside", whereas any non-zero pixel is
 considered outside. As a consequence, a pixel is considered to be on the boundary if its
 value is 0 and at least one of its 8-neighbour is not zero.
 
 "distance_transform" returns a pointer to an array of "width*height" floating point numbers
 containing for each pixel the distance to the closest boundary point.
 This array is allocated with `DISTANCE_TRANSFORM_MALLOC`. Use `DISTANCE_TRANSFORM_FREE`
 to delete it.
 
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

float* distance_transform(const uint8_t *in, int width, int height);

#ifdef __cplusplus
}
#endif

#endif /* DISTANCE_TRANSFORM_INCLUDE_H */

#ifdef DISTANCE_TRANSFORM_IMPLEMENTATION

#if defined(DISTANCE_TRANSFORM_MALLOC) && defined(DISTANCE_TRANSFORM_FREE)
// ok
#elif !defined(DISTANCE_TRANSFORM_MALLOC) && !defined(DISTANCE_TRANSFORM_FREE)
// ok
#else
#error "Must define all or none of DISTANCE_TRANSFORM_MALLOC and DISTANCE_TRANSFORM_FREE"
#endif

#if !defined(DISTANCE_TRANSFORM_MALLOC)
#define DISTANCE_TRANSFORM_MALLOC(sz) malloc(sz)
#define DISTANCE_TRANSFORM_FREE(p)    free(p)
#endif

float* distance_transform(const uint8_t *in, int width, int height) {
	int i, j;
	uint8_t c;
	const float max_dist = width*width + height*height;
	float *out = (float*)DISTANCE_TRANSFORM_MALLOC(width * height * sizeof(float));
	float *ptr = out;
	float *line;

	if(out == NULL) {
		return NULL;
	}

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
	ptr[0] = fmin(ptr[0], (fabs((dx)-(dy)) >= 1) ? (fmin((dx), (dy)) + 1.f) : (((dx) + (dy) + sqrt(2.f - ((dx)-(dy))*((dx)-(dy)))) / 2.f)); \
	ptr += (inc); \
} while(0)

	// y=[0,height[ x=[0,width[
	ptr = out;
	update_distance(ptr[1], ptr[width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[width], +1);
	}
	update_distance(ptr[-1], ptr[width], +1);

	line = ptr;
	for(j=1; j<(height-1); j++, line+=width) {
		ptr = line;
		update_distance(ptr[1], fmin(ptr[-width], ptr[width]), +1);
		for(i=1; i<(width-1); i++) {
			update_distance(fmin(ptr[-1], ptr[1]), fmin(ptr[-width], ptr[width]), +1);
		}
		update_distance(ptr[-1], fmin(ptr[-width], ptr[width]), +1);
	}

	update_distance(ptr[1], ptr[-width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[-width], +1);
	}
	update_distance(ptr[-1], ptr[-width], +1);

	// y=[0,height[ x=]width,0]
	ptr = out+width-1;
	update_distance(ptr[-1], ptr[width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[width], -1);
	}
	update_distance(ptr[1], ptr[width], -1);

	line = out+width;
	for(j=1; j<(height-1); j++, line+=width) {
        ptr = line+width-1;

		update_distance(ptr[-1], fmin(ptr[-width], ptr[width]), -1);
		for(i=1; i<(width-1); i++) {
			update_distance(fmin(ptr[-1], ptr[1]), fmin(ptr[-width], ptr[width]), -1);
		}
		update_distance(ptr[1], fmin(ptr[-width], ptr[width]), -1);
	}

	ptr = line + width-1;
	update_distance(ptr[-1], ptr[-width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[-width], -1);
	}
	update_distance(ptr[1], ptr[-width], -1);

	// y=]height,0] x=]width,0]
	ptr = out+width-1+(height-1)*width;
	update_distance(ptr[-1], ptr[-width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[-width], -1);
	}
	update_distance(ptr[1], ptr[-width], -1);

	for(j=1; j<(height-1); j++) {
		update_distance(ptr[-1], fmin(ptr[-width], ptr[width]), -1);
		for(i=1; i<(width-1); i++) {
			update_distance(fmin(ptr[-1], ptr[1]), fmin(ptr[-width], ptr[width]), -1);
		}
		update_distance(ptr[1], fmin(ptr[-width], ptr[width]), -1);
	}

	update_distance(ptr[-1], ptr[width], -1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[width], -1);
	}
	update_distance(ptr[1], ptr[width], -1);

	// y=]height,0] x=[0,width[
	line = ptr = out+(height-1)*width;
	update_distance(ptr[1], ptr[-width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[-width], +1);
	}
	update_distance(ptr[-1], ptr[-width], +1);

	line -= width;
	for(j=1; j<(height-1); j++, line-=width) {
		ptr = line;
		update_distance(ptr[1], fmin(ptr[-width], ptr[width]), +1);
		for(i=1; i<(width-1); i++) {
			update_distance(fmin(ptr[-1], ptr[1]), fmin(ptr[-width], ptr[width]), +1);
		}
		update_distance(ptr[-1], fmin(ptr[-width], ptr[width]), +1);
	}

	ptr = line;
	update_distance(ptr[1], ptr[width], +1);
	for(i=1; i<(width-1); i++) {
		update_distance(fmin(ptr[-1], ptr[1]), ptr[width], +1);
	}
	update_distance(ptr[-1], ptr[width], +1);

#undef update_distance

	return out;
}

#endif /* DISTANCE_TRANSFORM_IMPLEMENTATION */

#include <stdlib.h>

#include "stb_image.h"
#include "stb_image_write.h"

#define DISTANCE_TRANSFORM_IMPLEMENTATION
#include "../distance.h"

// Convert floating point value to RGB using hsl to rgb with s = 1 and l = 0.5
void heatmap(float s, float *r, float *g, float *b) {
	s = 1.0 - s;
	*r = fabs(fmod(s*6.0    , 6.0) - 3.0) - 1.0;
	*g = fabs(fmod(s*6.0+4.0, 6.0) - 3.0) - 1.0;
	*b = fabs(fmod(s*6.0+2.0, 6.0) - 3.0) - 1.0;
	*r = fmin(fmax(*r, 0.f), 1.f);
	*g = fmin(fmax(*g, 0.f), 1.f);
	*b = fmin(fmax(*b, 0.f), 1.f);
}

int main(int argc, char **argv) {
	uint8_t *input;
	int width, height;

	uint8_t *output;
	float *distance;

	if(argc != 3) {
		fprintf(stderr, "Usage: compute <input> <output>\n input - 8bpp greyscale image\n output - output istance map\n");
		return EXIT_FAILURE;
	}

	input = stbi_load(argv[1], &width, &height, NULL, 1);
	if(input == NULL) {
		fprintf(stderr, "Failed to read %s\n", argv[1]);
		return EXIT_FAILURE;
	}

	distance = (float*)malloc(width * height * sizeof(float));
	if(distance == NULL) {
		fprintf(stderr, "Failed to allocate distance map\n");
		return EXIT_FAILURE;
	}

	distance_transform(input, distance, width, height);
	
	// find the maximum distance in order to scale the values to [0,1]
	float max_dist = -1.f;
	float *ptr = distance;
	for(int j=0; j<height; j++) {
		for(int i=0; i<width; i++) {
			max_dist = fmax(max_dist, *ptr++);
		}
	}

	// generate "heatmap"
	output = (uint8_t*)malloc(width * height * 3);

	uint8_t *out = output;
	ptr = distance;
	for(int j=0; j<height; j++) {
		for(int i=0; i<width; i++) {
			float r, g, b, d = *ptr++ / max_dist;
			heatmap(d, &r, &g, &b);
			*out++ = r*255.f;
			*out++ = g*255.f;
			*out++ = b*255.f;
		}
	}

	stbi_write_png(argv[2], width, height, 3, output, 0);

	stbi_image_free(input);

	free(distance);

	free(output);

	return EXIT_SUCCESS;
}

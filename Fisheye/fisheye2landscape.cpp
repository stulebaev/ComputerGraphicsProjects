// fisheye2landscape.cpp
// cl /D_USE_MATH_DEFINES fisheye2landscape.cpp /MD jpeg.lib

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "jpeglib.h"


JSAMPLE* loadJPEG(const char* fileName, int& width, int& height) {
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  FILE* infile = fopen(fileName, "rb");
  if (infile == NULL) {
    fprintf(stderr, "Error opening jpeg file %s\n!", fileName);
    return NULL;
  }
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);
  jpeg_stdio_src(&cinfo, infile);
  jpeg_read_header(&cinfo, TRUE);
  width = cinfo.image_width;
  height = cinfo.image_height;
  jpeg_start_decompress(&cinfo);
  int row_stride = cinfo.output_width * cinfo.output_components;
  JSAMPARRAY buffer = (JSAMPARRAY)malloc(sizeof(JSAMPROW));
  buffer[0] = (JSAMPROW)malloc(sizeof(JSAMPLE)*row_stride);
  size_t image_size = row_stride * cinfo.output_height;
  JSAMPLE* image_data = (JSAMPLE*)malloc(image_size);
  unsigned long location = 0;
  while (cinfo.output_scanline < cinfo.output_height) {
    jpeg_read_scanlines(&cinfo, buffer, 1);
    memcpy(image_data+location, buffer[0], row_stride);
    location += row_stride;
  }
  jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);
  free(buffer[0]);
  fclose(infile);
  return image_data;
}

int saveJPEG(const char* fileName, int width, int height, unsigned char* imageData, int quality) {
  FILE* fout = fopen(fileName, "wb");
  if (fout == NULL) {
    fprintf(stderr, "Can't writing jpeg file %s\n!", fileName);
    return 0;
  }
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);
  jpeg_stdio_dest(&cinfo, fout);
  cinfo.image_width = width;
  cinfo.image_height = height;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;
  jpeg_set_defaults(&cinfo);
  jpeg_set_quality(&cinfo, quality, TRUE);
  jpeg_start_compress(&cinfo, TRUE);
  int row_stride = width*3;
  JSAMPROW row_pointer[1];
  while (cinfo.next_scanline < cinfo.image_height) {
    row_pointer[0] = &imageData[cinfo.next_scanline*row_stride];
    jpeg_write_scanlines(&cinfo, row_pointer, 1);
  }
  jpeg_finish_compress(&cinfo);
  fclose(fout);
  jpeg_destroy_compress(&cinfo);
  return 1;
}

bool fisheye2landscape(unsigned char* fisheye, int width, int height,
  unsigned char* landscape)
{
  /* the source image must be square and its width is even number */
  if ((width != height) || (width%2)) return false;

  int l = height/2;
  int w = 4*l;
  for (int i = 0; i < l; ++i) {
    for (int j = 0; j < w; ++j) {
      double radius = (double)(l - i);
      double theta = -2.0*M_PI*j/w;
      double fTrueX = radius*cos(theta);
      double fTrueY = radius*sin(theta);
      int x = (int)(fTrueX+0.5) + l;
      int y = l - (int)(fTrueY+0.5);
      /* check bounds */
      if (x >= 0 && x < width && y >= 0 && y < height) {
        /* getPixel(fisheye, x, y) */
        unsigned long location = (y*width + x)*3;
        unsigned char red = fisheye[location];
        location++;
        unsigned char green = fisheye[location];
        location++; 
        unsigned char blue = fisheye[location];
        /* setPixel(landscape, j, i) */
        location = (i*w + j)*3;
        landscape[location] = red;
        location++;
        landscape[location] = green;
        location++;
        landscape[location] = blue;
      }
    }
  }
  return true;
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: fisheye2landscape input_file\n");
    return -1;
  }

  int width, height;
  unsigned char* fisheye;
  if ((fisheye = loadJPEG(argv[1], width, height)) == NULL)
    return -1;

  unsigned char* landscape;
  landscape = (unsigned char*)malloc(width*width*3);
  if (landscape == NULL) {
    free(fisheye);
    return -1;
  }

  if (fisheye2landscape(fisheye, width, height, landscape)) {
    char outName[FILENAME_MAX];
    char* ext;
    strcpy(outName, argv[1]);
    ext = strrchr(outName, '.');
    strcpy(ext, "_landscape.jpg");
    saveJPEG(outName, 2*width, height/2, landscape, 75);
  }
  free(fisheye);
  free(landscape);

  return 0;
}

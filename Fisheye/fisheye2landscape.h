#ifndef FISHEYE2LANDSCAPE_H
#define FISHEYE2LANDSCAPE_H

class VideoProcess {
public:
  VideoProcess(void* imagePointer, int width, int height);
  ~VideoProcess();
  bool resultInfo(void* buffer, int& targetWidth, int& targetHeight) const;
  bool process();
};

#endif //FISHEYE2LANDSCAPE_H
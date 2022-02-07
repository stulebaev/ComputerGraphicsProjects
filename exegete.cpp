// exegete.cpp

#include <list>
#include <vector>
#include <string>
#include <iostream>
#include <Magick++.h>

std::list<Magick::Image> imagesList;

void tokenize(std::string& str, const char delim, std::vector<std::string>& out)
{
  size_t start, end = 0;
  while ((start = str.find_first_not_of(delim, end)) != std::string::npos) {
    end = str.find(delim, start);
    out.push_back(str.substr(start, end - start));
  }
}

void print_help(void)
{
  std::cout << "Supported commands:" << std::endl;
  std::cout << '\t' << "load <name> <filename>" << std::endl;
  std::cout << '\t' << "save <name> <filename>" << std::endl;
  std::cout << '\t' << "blur <from_name> <to_name> <size>" << std::endl;
  std::cout << '\t' << "resize <from_name> <to_name> <new_width> <new_height>" << std::endl;
}

void load_image(std::string& name, std::string& filename)
{
  try {
    Magick::Image image(filename);
    image.label(name);
    imagesList.push_back(image);
    std::cout << "Image '" << name << "' was loaded from file '" << filename << "'" << std::endl;
  }
  catch (Magick::Exception& error) { 
    std::cerr << error.what() << std::endl;
  }
}

void save_image(std::string& name, std::string& filename)
{
  for (auto& image : imagesList)
    if (image.label() == name) {
      try {
        image.write(filename);
        std::cout << "Image '" << name << "' was saved to file '" << filename << "'" << std::endl;
        return;
      }
      catch (Magick::Exception& error) { 
        std::cerr << error.what() << std::endl;
      }
    }
  std::cout << "Image '" << name << "' not found" << std::endl;
}

void blur_image(std::string& from_name, std::string& to_name, std::string& size)
{
  double radius;
  try { radius = std::stod(size); } catch (...) {
    std::cerr << "Invalid size" << std::endl;
    return;
  }
  for (auto& image : imagesList)
    if (image.label() == from_name) {
      try {
        Magick::Image blurred(image);
        blurred.blur(radius);
        blurred.label(to_name);
        imagesList.push_back(blurred);
        std::cout << "Image '" << from_name << "' was blurred to image '" << to_name << "'" << std::endl;
        return;
      }
      catch (Magick::Exception& error) { 
        std::cerr << error.what() << std::endl;
      }
    }
  std::cout << "Image '" << from_name << "' not found" << std::endl;
}

void resize_image(std::string& from_name, std::string& to_name,
                  std::string& new_width, std::string& new_height)
{
  size_t width_, height_;
  try { width_ = std::stoul(new_width); } catch (...) {
    std::cerr << "Invalid new width" << std::endl;
    return;
  }
  try { height_ = std::stoul(new_height); } catch (...) {
    std::cerr << "Invalid new height" << std::endl;
    return;
  }
  for (auto& image : imagesList)
    if (image.label() == from_name) {
      try {
        Magick::Image resized(image);
        resized.resize(Magick::Geometry(width_,height_));
        resized.label(to_name);
        imagesList.push_back(resized);
        std::cout << "Image '" << from_name << "' was resized to image '" << to_name << "'" << std::endl;
        return;
      }
      catch (Magick::Exception& error) { 
        std::cerr << error.what() << std::endl;
      }
    }
  std::cout << "Image '" << from_name << "' not found" << std::endl;
}

int main(int argc, char *argv[])
{
  const char* PROMPT = ">";

  Magick::InitializeMagick(*argv);

  while (true) {
    std::string line;
    std::vector<std::string> tokens;

    std::cout << PROMPT;
    std::getline(std::cin, line);
    tokenize(line, ' ', tokens);

    if ( tokens[0] == "quit" ||
         tokens[0] == "exit" ||
         tokens[0].front() == 'q' ) break;

    if ( tokens[0] == "help" ||
         tokens[0].front() == 'h' )
      print_help();

    if ( tokens[0] == "load" ||
         tokens[0].front() == 'l' )
      if (tokens.size() == 3) load_image(tokens[1], tokens[2]);

    if ( tokens[0] == "save" ||
         tokens[0].front() == 's' )
      if (tokens.size() == 3) save_image(tokens[1], tokens[2]);

    if ( tokens[0] == "blur" )
      if (tokens.size() == 4) blur_image(tokens[1], tokens[2], tokens[3]);

    if ( tokens[0] == "resize" )
      if (tokens.size() == 5) resize_image(tokens[1], tokens[2], tokens[3], tokens[4]);
  }

  return 0;
}

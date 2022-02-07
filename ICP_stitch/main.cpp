#include <vector>
#include <string>
#include <cstdio>
#include <boost/filesystem.hpp>
#include "ply_io.h"
#include "icp.h"

PlyProperty vert_props[] = { /* list of property information for a vertex */
  {"x", Float32, Float32, offsetof(Point3f,X), 0, 0, 0, 0},
  {"y", Float32, Float32, offsetof(Point3f,Y), 0, 0, 0, 0},
  {"z", Float32, Float32, offsetof(Point3f,Z), 0, 0, 0, 0},
};

void savePLY(std::string filename, std::vector<Point3f>& vertices)
{
	unsigned int numVertices = vertices.size();

	// Open File
	FILE *meshFile = fopen(filename.c_str(), "wt");

	// Write the header line
	std::string header = "ply\nformat ascii 1.0\n";
	fwrite(header.c_str(), sizeof(char), header.length(), meshFile);

	const unsigned int bufSize = 1000 * 3;
	char outStr[bufSize];
	int written = 0;

	// Elements are: x,y,z
	written = sprintf(outStr, "element vertex %u\nproperty float x\nproperty float y\nproperty float z\n", numVertices);
	fwrite(outStr, sizeof(char), written, meshFile);

	written = sprintf(outStr, "end_header\n");
	fwrite(outStr, sizeof(char), written, meshFile);

	for (unsigned int vertexIndex = 0; vertexIndex < numVertices; vertexIndex++)
	{
		written = sprintf(outStr, "%f %f %f\n", vertices[vertexIndex].X, vertices[vertexIndex].Y, vertices[vertexIndex].Z);
		fwrite(outStr, sizeof(char), written, meshFile);
	}

	fflush(meshFile);
	fclose(meshFile);
}

int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		fprintf(stderr, "Usage: %s  directory_with_PLY_files\n", argv[0]);
		return -1;
	}

	boost::filesystem::path directory(argv[1]);
	if (!boost::filesystem::is_directory(directory))
	{
		fprintf(stderr, "Error: %s is not a directory\n", directory.string().c_str());
		return -1;
	}

	std::vector<std::string> matching_files;
	boost::filesystem::directory_iterator end_itr;
	for (boost::filesystem::directory_iterator itr(directory); itr != end_itr; ++itr)
	{
		// Skip if not a file
		if (!boost::filesystem::is_regular_file(itr->status())) continue;
		if (itr->path().extension() == ".ply")
			// File matches, store it
			matching_files.push_back(itr->path().string());
	}
	if (matching_files.size() == 0)
	{
		fprintf(stderr, "Error: PLY files are not found in %s\n", directory.string().c_str());
		return -1;
	}

	std::vector<Point3f> vertices1, vertices2;
	for (auto& filename : matching_files)
	{
		//printf("...reading file %s\n", filename.c_str());
		PlyFile* plyfile;
		FILE* fp = fopen(filename.c_str(), "r");
		if ((fp == NULL) || ((plyfile = read_ply(fp)) == NULL))
		{
			fprintf(stderr, "Error loading point cloud from file '%s'\n", filename.c_str());
			break;
		}
		int elem_count;
		char* elem_name;
		for (int i = 0; i < plyfile->num_elem_types; i++)
		{
			elem_name = setup_element_read_ply(plyfile, i, &elem_count);
			if (equal_strings("vertex", elem_name))
			{
				setup_property_ply(plyfile, &vert_props[0]);
				setup_property_ply(plyfile, &vert_props[1]);
				setup_property_ply(plyfile, &vert_props[2]);
				for (int j = 0; j < elem_count; j++)
				{
					Point3f vertex;
					get_element_ply(plyfile, &vertex);
					vertices1.push_back(vertex);
				}
			}
		}
		fclose(plyfile->fp);
		free(plyfile);

		if (vertices2.empty())
		{
			vertices2 = vertices1;
			continue;
		}
		cv::Mat R = cv::Mat::eye(3, 3, CV_32F);
		cv::Mat t(1, 3, CV_32F, cv::Scalar(0));
		ICP(vertices1.data(), vertices1.size(), vertices2.data(), vertices2.size(), (float*)R.data, (float*)t.data, 1);
		vertices1.clear();
	}
	//for (auto& vertex : vertices2) printf("%f %f %f\n", vertex.X, vertex.Y, vertex.Z);
	savePLY("./testResult.ply", vertices2);

	return 0;
}

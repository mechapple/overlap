/*
 * read_STL.cpp
 * 
 * Copyright 2017 Anirban <anirban@ZeroPointEnergy>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 * 
 */


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <stdlib.h>
#include <map>

#define TOL2 1e-12

struct triplet {
	int v1,v2,v3;
	triplet(int a, int b, int c) : v1(a), v2(b), v3(c) {}
};

struct point2 {
  double x,y,z;
  point2(double a, double b, double c) : x(a), y(b), z(c) {}
  
  bool operator <(const point2& pt) const
  {
	return (x < pt.x) || ((!(pt.x < x)) && (y < pt.y)) || ((!(pt.x < x)) && (!(pt.y < y)) && (z < pt.z));
  }
};

int main(int argc, char **argv)
{
	std::ifstream infile0(argv[1]);
	std::string line,input;
	std::stringstream sstream;
	int count=0,vertex_count=0,face_count=0;;
	
	std::map<point2, int> vertex_points;
	
	while(std::getline(infile0, line))
	{
		std::vector<std::string> entries;
		sstream.clear(); sstream.str(line);
		while(sstream >> input)
		{
			entries.push_back(input);
		}
		if(entries[0]=="vertex") 
		{
			point2 p1(::atof(entries[1].c_str()),::atof(entries[2].c_str()),::atof(entries[3].c_str()));
			
			if(vertex_points.find(p1) == vertex_points.end())
			{
				vertex_points.insert(std::pair<point2,int>(p1,vertex_count++));
			}	
			count++;
		}
	}
	infile0.close();
	
	face_count=count/3;
	//std::cout << "\n Numvertices = " << vertex_count << "\n";
	//std::cout << "\n Numfaces = " << face_count << "\n";
	std::vector<point2> vertices; for(int i=0;i<vertex_count;i++) vertices.push_back(point2(0,0,0));
	std::vector<triplet> faces;
	
	for(std::map<point2, int>::const_iterator it = vertex_points.begin(); it != vertex_points.end(); ++it)
	{
		//int vloc = it->second;
		vertices[it->second].x = it->first.x; vertices[it->second].y = it->first.y; vertices[it->second].z = it->first.z;
		//std::cout << it->first << " " << it->second << " " << Svertices[it->second].x << " " << Svertices[it->second].y << " " << Svertices[it->second].z << "\n";
	}
	
	count=0;
	std::ifstream infile1(argv[1]);
	while(std::getline(infile1, line))
	{
		std::vector<std::string> entries;
		sstream.clear(); sstream.str(line);
		while(sstream >> input)
		{
			entries.push_back(input);
		}
		
		if(entries[0]=="outer") 
		{
			int vloc[3];
			for(int i=0;i<3;i++)
			{
				std::getline(infile1, line);
				
				std::vector<std::string> entries1;
				sstream.clear(); sstream.str(line);
				while(sstream >> input)
				{
					entries1.push_back(input);
				}
			
				point2 p1(::atof(entries1[1].c_str()),::atof(entries1[2].c_str()),::atof(entries1[3].c_str()));
				vloc[i] = vertex_points.at(p1);
			}
			
			faces.push_back(triplet(vloc[0],vloc[1],vloc[2]));
			//point2 p1(::atof(entries[1].c_str()),::atof(entries[2].c_str()),::atof(entries[3].c_str()));
			
			//if(vertex_points.find(p1) == vertex_points.end())
			//{
				//vertex_points.insert(std::pair<point2,int>(p1,vertex_count++));
			//}	
			count++;
		}
	}
	
	std::cout << "\n Numvertices = " << vertices.size() << "\n";
	std::cout << "\n Numfaces = " << faces.size() << "\n";
	
	infile1.close();
	
	std::ofstream outfileV("Vertices.dat");
	std::ofstream outfileF("Faces.dat");
	
	for(int i=0;i<vertices.size();i++)
		outfileV << vertices[i].x << " " << vertices[i].y << " " << vertices[i].z << "\n";
	
	for(int i=0;i<faces.size();i++)
		outfileF << faces[i].v1 << " " << faces[i].v2 << " " << faces[i].v3 << "\n";
	
	outfileV.close();
	outfileF.close();
	
	return 0;
}


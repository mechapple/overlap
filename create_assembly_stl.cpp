/*  g++ -L/home/anirban/Downloads/boost_1_59_0 -L/home/anirban/Downloads/boost_1_59_0/stage/lib 
 * 				-g -O3 -o exec.o file.cpp -lboost_timer -lboost_system
 * create_Si_system.cpp
 * 
 * THIS creates an LMP file and STL file. 
 * USAGE: ./create_assembly_stl_lmp <fraction> <size> <percolated yes/no(1/0)>
 * 
 * THIS also adds internal faces to the components but not their inverse
 * 
 * Copyright 2015 Anirban <anirban@ZeroPointEnergy>
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
 * aint with this program; if not, write to the Free Software
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
#include <stack>
#include <cmath>
#include <stdlib.h>
#include <set>
#include <map>
#include <ctime>

#include <boost/config.hpp>
#include <algorithm>
#include <utility>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/graph_traits.hpp>

#include <boost/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/timer/timer.hpp>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/edmonds_karp_max_flow.hpp>

#include <boost/filesystem.hpp>

#define TOL 0.00 //TOL must be ZERO for this to work well
#define RTOL 0.001
#define EPS 1e-16
#define EPS2 1e-12

#define TOL2 1e-6
#define lata 4.0000

typedef double_t EdgeWeightType;


typedef boost::adjacency_list_traits < boost::vecS, boost::vecS, boost::directedS > Traits;

typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::directedS,
  boost::property < boost::vertex_name_t, std::string,
    boost::property < boost::vertex_index_t, int,
      boost::property < boost::vertex_color_t, boost::default_color_type,
        boost::property < boost::vertex_distance_t, double_t,
          boost::property < boost::vertex_predecessor_t, Traits::edge_descriptor > > > > >,
 
  boost::property < boost::edge_capacity_t, EdgeWeightType,
    boost::property < boost::edge_residual_capacity_t, EdgeWeightType,
      boost::property < boost::edge_reverse_t, Traits::edge_descriptor > > > > Graph;

Traits::edge_descriptor AddEdge(Traits::vertex_descriptor &v1,
                                Traits::vertex_descriptor &v2,
                                boost::property_map < Graph, boost::edge_reverse_t >::type &rev,
                                const double capacity,
                                Graph &g);
                                
Traits::edge_descriptor AddEdge(Traits::vertex_descriptor &v1, Traits::vertex_descriptor &v2, boost::property_map < Graph, boost::edge_reverse_t >::type &rev, const double capacity, Graph &g)
{
	Traits::edge_descriptor e1 = add_edge(v1, v2, g).first;
	Traits::edge_descriptor e2 = add_edge(v2, v1, g).first;
	put(boost::edge_capacity, g, e1, capacity);
	put(boost::edge_capacity, g, e2, capacity);
	
	rev[e1] = e2;
	rev[e2] = e1;
}

struct point {
  double x,y,z;
  int id,color;
  int mol_id;
  std::vector<int> sfaces,bfaces,edges;
  std::set<int> grains;
  double pbc_disp[3];
  point(double a, double b, double c, int d, int e) : x(a), y(b), z(c), id(d), mol_id(e) {}
};

struct point2 {
  double x,y,z;
  point2(double a, double b, double c) : x(a), y(b), z(c) {}
  
  bool operator <(const point2& pt) const
  {
	return (x < pt.x) || ((!(pt.x < x)) && (y < pt.y)) || ((!(pt.x < x)) && (!(pt.y < y)) && (z < pt.z));
  }
};

struct edge {
  //double x,y,z;
  point *m,*n;
  int id,p,q,active,pbc_face,pbc_zface;
  std::vector<int> plane_face;
  double len;
  double normal[3],normal2[3],center[3];
  std::vector<double> vertices;
  std::vector<std::string> vertex_strs;
  edge(point *a, point *b, int d, double l) : m(a), n(b), id(d), len(l) 
  {
	pbc_face=0; 
	pbc_zface=0;
	//plane_face=-1;
  }
};

struct complist {
	std::vector<int> contents;
};

struct boundary {
	double xlo,xhi;
	double ylo,yhi;
	double zlo,zhi;
	double xy,xz,yz;
};

void minimum_image(double &dx, double &dy, double &dz, boundary b)
{
	//printf("\n %g %g %g %g %g %g %g %g %g\n",b.xlo,b.xhi,b.ylo,b.yhi,b.zlo,b.zhi,b.xy,b.xz,b.yz);
	double xprd=(b.xhi-b.xlo), yprd=(b.yhi-b.ylo), zprd=(b.zhi-b.zlo);
	
	if (fabs(dz) > zprd*0.5) {
		if (dz < 0.0) {
			dz += zprd;
			dy += b.yz;
			dx += b.xz;
		} else {
			dz -= zprd;
			dy -= b.yz;
			dx -= b.xz;
		}
	}

	if (fabs(dy) > yprd*0.5) {
		if (dy < 0.0) {
			dy += yprd;
			dx += b.xy;
		} else {
			dy -= yprd;
			dx -= b.xy;
		}
	}
	
	if (fabs(dx) > xprd*0.5) {
		if (dx < 0.0) dx += xprd;
		else dx -= xprd;
	}
}

void closest_image(double &dx, double &dy, double &dz, boundary b)
{
	//double dx = xj[0] - xi[0];
	//double dy = xj[1] - xi[1];
	//double dz = xj[2] - xi[2];
	double xprd=(b.xhi-b.xlo), yprd=(b.yhi-b.ylo), zprd=(b.zhi-b.zlo);
	
	if (dz < 0.0) {
		while (dz < 0.0) {
			dz += zprd;
			dy += b.yz;
			dx += b.xz;
		}
		if (dz > zprd*0.5) {
			dz -= zprd;
			dy -= b.yz;
			dx -= b.xz;
		}
	} else {
		while (dz > 0.0) {
			dz -= zprd;
			dy -= b.yz;
			dx -= b.xz;
		}
		if (dz < -zprd*0.5) {
			dz += zprd;
			dy += b.yz;
			dx += b.xz;
		}
	}
    

	if (dy < 0.0) {
		while (dy < 0.0) {
			dy += yprd;
			dx += b.xy;
		}
		if (dy > yprd*0.5) {
			dy -= yprd;
			dx -= b.xy;
		}
	} else {
		while (dy > 0.0) {
			dy -= yprd;
			dx -= b.xy;
		}
		if (dy < -yprd*0.5) {
			dy += yprd;
			dx += b.xy;
		}
	}
    
    if (dx < 0.0) {
		while (dx < 0.0) dx += xprd;
		if (dx > xprd*0.5) dx -= xprd;
	} else {
		while (dx > 0.0) dx -= xprd;
		if (dx < -xprd*0.5) dx += xprd;
	}
    
	//xjimage[0] = xi[0] + dx;
	//xjimage[1] = xi[1] + dy;
	//xjimage[2] = xi[2] + dz;
}
/*
double rn(boost::mt19937 rng)
{
  //boost::mt19937 rng(time(0));
  //static boost::uniform_01<boost::mt19937> zeroone(rng);
  
  return zeroone();
}*/

int main(int argc, char **argv)
{
	int start_s=clock();
	
	//int COMP_CUTOFF = 3;
	//if(argc>3) COMP_CUTOFF = atoi(argv[3]);
	if(argc<3) {
		printf("\n USAGE: ./create_assembly_stl <num_mols> <nx> <ny> <nz>"); 
		exit(EXIT_FAILURE);
	}
	
	double direction[3] = {0.0,0.0,1.0};
		
	std::string line,thrash,thrash1;
	int atomspmol = atoi(argv[1]);
	int nx = atoi(argv[2]),ny = atoi(argv[3]),nz = atoi(argv[4]);
	
	boundary bounds;
	
	//const int RANDOM_SEED = 2806;
	//boost::mt19937 rng(time(0));
	long int RANDOM_SEED = time(0);
	std::cout << "\nRandom Seed = " << RANDOM_SEED;
	boost::mt19937 rng(RANDOM_SEED);
	static boost::uniform_01<boost::mt19937> rn(rng);
	
	//##// 1. Read Voronoi cell data from file
	
	std::ifstream infile0("voro_cell.dat");
	int count=0;
	while(std::getline(infile0, line)) count++;
	int ngrains = count, nedges = 0, nbedges = 0;
	//std::cout << count;
		
	std::stringstream sstream;
	
	std::vector<point> Sgrains,COMS;
	std::vector<edge> Sedges,Bedges;
	std::vector<int> leftF,rightF;
	
	std::set<double> plane_set;
	
	for(int i=0;i<ngrains;i++) Sgrains.push_back(point(0.00001,0.00001,0.00001,i+1,(i/atomspmol)));
	
	//for(int i=0;i<ngrains;i++) std::cout << std::endl << Sgrains[i].id << " " << Sgrains[i].x << " " << Sgrains[i].y << " " << Sgrains[i].z;
			
	std::string input;
	
	infile0.close();
	
	//******** BOUNDARY INFO ***********//
	std::ifstream infileb("bounds.dat");
	
	//double xlo=0.0,xhi=1.0;
	//double ylo=0.0,yhi=1.0;
	//double zlo=0.0,zhi=1.0;
	//double xy=0.0,xz=0.0,yz=0.0;
	int countb=0,z_switch=0;
	
	while(std::getline(infileb, line))
	{
		//outLfile << line << std::endl;
		std::vector<std::string> entries;
		sstream.clear(); sstream.str(line);
		while(sstream >> input)
		{
			entries.push_back(input);
		}
		countb++;
		if(countb==1) { bounds.xlo = ::atof(entries[0].c_str()); bounds.xhi = ::atof(entries[1].c_str());}
		if(countb==2) { bounds.ylo = ::atof(entries[0].c_str()); bounds.yhi = ::atof(entries[1].c_str());}
		if(countb==3) { bounds.zlo = ::atof(entries[0].c_str()); bounds.zhi = ::atof(entries[1].c_str());}
		if(countb==4) { bounds.xy = ::atof(entries[0].c_str()); bounds.xz = ::atof(entries[1].c_str()); bounds.yz = ::atof(entries[2].c_str());}
	}
	
	printf("\n %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g\n",bounds.xlo,bounds.xhi,bounds.ylo,bounds.yhi,bounds.zlo,bounds.zhi,bounds.xy,bounds.xz,bounds.yz);
	
	//double xprd=(xhi-xlo), yprd=(yhi-ylo), zprd=(zhi-zlo);
	infileb.close();
	//**********************************//
	
	//******** COMS INFO ***********//
	std::ifstream infilec("hmx_coms.dat");
	while(std::getline(infilec, line))
	{
		std::vector<std::string> entries;
		sstream.clear(); sstream.str(line);
		while(sstream >> input)
		{
			entries.push_back(input);
		}
		
		int mol_id = atoi(entries[0].c_str())-1;
		COMS.push_back(point(::atof(entries[1].c_str()),::atof(entries[2].c_str()),::atof(entries[3].c_str()),mol_id,1));
		plane_set.insert(::atof(entries[3].c_str()));
	}
	
	if(COMS[0].z<COMS[1].z) z_switch=0;
	else z_switch=1;
	
	infilec.close();
	//**********************************//
	
	std::vector<double> COMplanes(plane_set.begin(), plane_set.end());
	std::vector<double> crackplanes(COMplanes.size()-1),crackplanes_zmin(COMplanes.size()-1);
	std::vector< std::vector<int> > crackplane_faces(COMplanes.size()-1);
	
	//int plane_id=0;
	//for(std::set<double>::const_iterator it = plane_set.begin(); it != plane_set.end(); ++it,plane_id++)
	//{
		//printf("\n Plane %d: %4.10g",plane_id+1,COMplanes[plane_id]);
		//printf("\n Plane %d: %4.10g",plane_id+1,*it);
	//}
	
	for(int i=0;i<crackplanes.size();i++)
	{
		crackplanes[i] = 0.5*(COMplanes[i]+COMplanes[i+1]);
		printf("\n crack Plane %d: %4.10g",i+1,crackplanes[i]);
	}
	
	std::ifstream infile1("voro_cell.dat");
	while(std::getline(infile1, line))
	{
		std::vector<std::string> entries;
		sstream.clear(); sstream.str(line);
		while(sstream >> input)
		{
			entries.push_back(input);
		}
		
		
		//for(int i=0; i< entries.size(); i++) std::cout << entries[i] << " "; std::cout << std::endl;
		
		int id_curr = atoi(entries[0].c_str())-1; //radical voronoi indexes from 1 to n
		Sgrains[id_curr].id = id_curr;
		Sgrains[id_curr].x = ::atof(entries[1].c_str()); Sgrains[id_curr].y = ::atof(entries[2].c_str()); Sgrains[id_curr].z = ::atof(entries[3].c_str());
		
		double dx0=Sgrains[id_curr].x-COMS[Sgrains[id_curr].mol_id].x, dy0=Sgrains[id_curr].y-COMS[Sgrains[id_curr].mol_id].y, dz0=Sgrains[id_curr].z-COMS[Sgrains[id_curr].mol_id].z;
		double dx=Sgrains[id_curr].x-COMS[Sgrains[id_curr].mol_id].x, dy=Sgrains[id_curr].y-COMS[Sgrains[id_curr].mol_id].y, dz=Sgrains[id_curr].z-COMS[Sgrains[id_curr].mol_id].z;
		
		closest_image(dx,dy,dz,bounds);
		
		//Sgrains[id_curr].pbc_disp[0] = 0.0; Sgrains[id_curr].pbc_disp[1] = 0.0; Sgrains[id_curr].pbc_disp[2] = 0.0;
		Sgrains[id_curr].pbc_disp[0] = (dx-dx0)*(fabs(dx-dx0)>EPS2); Sgrains[id_curr].pbc_disp[1] = (dy-dy0)*(fabs(dy-dy0)>EPS2); Sgrains[id_curr].pbc_disp[2] = (dz-dz0)*(fabs(dz-dz0)>EPS2);
		
		Sgrains[id_curr].x += Sgrains[id_curr].pbc_disp[0];
		Sgrains[id_curr].y += Sgrains[id_curr].pbc_disp[1];
		Sgrains[id_curr].z += Sgrains[id_curr].pbc_disp[2];
		
		//std::cout << std::endl << id_curr << " " << Sgrains[id_curr].id << " " << Sgrains[id_curr].pbc_disp[0] << " " << Sgrains[id_curr].pbc_disp[1] << " " << Sgrains[id_curr].pbc_disp[2];
		//printf("\n %d %d %4.16g %4.16g %4.16g",id_curr,Sgrains[id_curr].id,Sgrains[id_curr].pbc_disp[0],Sgrains[id_curr].pbc_disp[1],Sgrains[id_curr].pbc_disp[2]);
		
		//std::cout << std::endl << id_curr << " " << Sgrains[id_curr].id << " " << Sgrains[id_curr].x << " " << Sgrains[id_curr].y << " " << Sgrains[id_curr].z;
		
		//COMS[Sgrains[id_curr].mol_id;
		
	}
	
	infile1.close();
	
	//*****************************
	std::ifstream infile("voro_cell.dat");
	
	while(std::getline(infile, line))
	{
		//outLfile << line << std::endl;
		std::vector<std::string> entries;
		sstream.clear(); sstream.str(line);
		while(sstream >> input)
		{
			entries.push_back(input);
		}
		
		//for(int i=0; i< entries.size(); i++) std::cout << entries[i] << " "; std::cout << std::endl;
		
		int id_curr = atoi(entries[0].c_str())-1; //radical voronoi indexes from 1 to n
		//Sgrains[id_curr].id = id_curr;
		//Sgrains[id_curr].x = ::atof(entries[1].c_str()); Sgrains[id_curr].y = ::atof(entries[2].c_str()); Sgrains[id_curr].z = ::atof(entries[3].c_str());
		
		//std::cout << std::endl << id_curr << " " << Sgrains[id_curr].id << " " << Sgrains[id_curr].x << " " << Sgrains[id_curr].y << " " << Sgrains[id_curr].z;
		
		
		int numfaces = atoi(entries[4].c_str());
		for(int i=5;i<(5+numfaces);i++)
		{
			int id2 = atoi(entries[i].c_str())-1;
			double face_area = ::atof(entries[i+numfaces].c_str());
			
			if (id2!=id_curr) 
			{
				Sedges.push_back(edge(&Sgrains[id_curr],&Sgrains[id2],nedges++,face_area));
				Sgrains[id_curr].sfaces.push_back(nedges-1); //only id_curr as per (id2!=id_curr) condition
				//Sgrains[id2].sfaces.push_back(nedges-1);
				
				std::string curr_string = entries[i+2*numfaces];
				curr_string = curr_string.substr(1,curr_string.size()-2);
				std::stringstream ss(curr_string);
				std::vector<int> face_vertex_ids;
				
				while( ss.good() )
				{
				    std::string substr;
				    std::getline(ss,substr,',');
				    face_vertex_ids.push_back(atoi(substr.c_str()));
				}
				
				for(int j=0;j<face_vertex_ids.size();j++)
				{
					//if(j<(face_vertex_ids.size()-1)) std::cout << "|";
					//else std::cout << " ";
					int local_vid = face_vertex_ids[j];
					curr_string = entries[local_vid+6+3*numfaces];
					curr_string = curr_string.substr(1,curr_string.size()-2);
					std::stringstream ss2(curr_string);
					
					Sedges.back().vertex_strs.push_back(entries[local_vid+6+3*numfaces]);
					//std::cout << "\n" << entries[local_vid+6+3*numfaces];
					int dim=0;
					while( ss2.good() )
					{
						std::string substr;
					    std::getline(ss2,substr,',');
					    
					    double vloc = ::atof(substr.c_str());
					    Sedges.back().vertices.push_back(vloc+Sgrains[id_curr].pbc_disp[dim]);
					    dim++;
					    //result.push_back(::atof(substr.c_str()));
					}
					
					//std::cout << face_vertex_ids[j] << " [" << curr_string << "] ";
					//if(j<(face_vertex_ids.size()-1)) std::cout << "|";
					//else std::cout << " ";
						
				}
				
				double face_center[3]={0.0,0.0,0.0};
				
				for(int j=0;j<face_vertex_ids.size();j++)
				{
					face_center[0] += Sedges.back().vertices[j*3];
					face_center[1] += Sedges.back().vertices[j*3+1];
					face_center[2] += Sedges.back().vertices[j*3+2];
				}
				
				//std::cout << "\n facevertex_num " << face_vertex_ids.size();
				//std::cout << "\n face_center0 " << face_center[0] << " " << face_center[1] << " " << face_center[2];
				
				for(int j=0;j<3;j++) face_center[j] = face_center[j]/face_vertex_ids.size();
				for(int j=0;j<3;j++) Sedges.back().center[j] = face_center[j];
				
				
				double dist1[3] = {face_center[0]-Sgrains[id_curr].x,face_center[1]-Sgrains[id_curr].y,face_center[2]-Sgrains[id_curr].z};
				double dist2[3] = {face_center[0]-Sgrains[id2].x,face_center[1]-Sgrains[id2].y,face_center[2]-Sgrains[id2].z};
				
				Sedges.back().normal[0] = dist1[0]; Sedges.back().normal[1] = dist1[1]; Sedges.back().normal[2] = dist1[2];
				
				double xprd=(bounds.xhi-bounds.xlo), yprd=(bounds.yhi-bounds.ylo), zprd=(bounds.zhi-bounds.zlo);
				
				//std::cout << "\n face_center " << face_center[0] << " " << face_center[1] << " " << face_center[2];
				//std::cout << "\n dist1 " << dist1[0] << " " << dist1[1] << " " << dist1[2];
				//std::cout << "\n dist2 " << dist2[0] << " " << dist2[1] << " " << dist2[2];
				
				int cond1 = (fabs(dist1[0])<xprd*0.5)&&(fabs(dist1[1])<yprd*0.5)&&(fabs(dist1[2])<zprd*0.5);
				int cond2 = (fabs(dist2[0])<xprd*0.5)&&(fabs(dist2[1])<yprd*0.5)&&(fabs(dist2[2])<zprd*0.5);
				
				//if(cond1==0 || cond2==0) Sedges.back().pbc_face=1;
				if(cond1==1 && cond2==0) Sedges.back().pbc_face=1;
				
				if(fabs(dist1[2])<zprd*0.5 && fabs(dist2[2])>zprd*0.5) Sedges.back().pbc_zface=1;
				
				//Sgrains[id_curr].sfaces.push_back(nedges-1);
			}
			
			if (id2==id_curr)  //when grain adjoins itself via PBC
			{
				Sedges.push_back(edge(&Sgrains[id_curr],&Sgrains[id2],nedges++,face_area));
				Sgrains[id_curr].sfaces.push_back(nedges-1);
				//Sgrains[id2].sfaces.push_back(nedges-1);
				
				Sedges.back().pbc_face=1;
				
				std::string curr_string = entries[i+2*numfaces];
				curr_string = curr_string.substr(1,curr_string.size()-2);
				std::stringstream ss(curr_string);
				std::vector<int> face_vertex_ids;
				
				int dim=0;
				while( ss.good() )
				{
				    std::string substr;
				    std::getline(ss,substr,',');
				    face_vertex_ids.push_back(atoi(substr.c_str()));
				}
				
				for(int j=0;j<face_vertex_ids.size();j++)
				{
					//if(j<(face_vertex_ids.size()-1)) std::cout << "|";
					//else std::cout << " ";
					int local_vid = face_vertex_ids[j];
					curr_string = entries[local_vid+6+3*numfaces];
					curr_string = curr_string.substr(1,curr_string.size()-2);
					std::stringstream ss2(curr_string);
					
					Sedges.back().vertex_strs.push_back(entries[local_vid+6+3*numfaces]);
					//std::cout << "\n" << entries[local_vid+6+3*numfaces];
					
					while( ss2.good() )
					{
						std::string substr;
					    std::getline(ss2,substr,',');
					    
					    double vloc = ::atof(substr.c_str());
					    Sedges.back().vertices.push_back(vloc+Sgrains[id_curr].pbc_disp[dim]);
					    dim++;
					}
					
					//std::cout << face_vertex_ids[j] << " [" << curr_string << "] ";
					//if(j<(face_vertex_ids.size()-1)) std::cout << "|";
					//else std::cout << " ";
						
				}
				
				double face_center[3]={0.0,0.0,0.0};
				
				for(int j=0;j<face_vertex_ids.size();j++)
				{
					face_center[0] += Sedges.back().vertices[j*3];
					face_center[1] += Sedges.back().vertices[j*3+1];
					face_center[2] += Sedges.back().vertices[j*3+2];
				}
				
				//std::cout << "\n facevertex_num " << face_vertex_ids.size();
				//std::cout << "\n face_center0 " << face_center[0] << " " << face_center[1] << " " << face_center[2];
				
				for(int j=0;j<3;j++) face_center[j] = face_center[j]/face_vertex_ids.size();
				for(int j=0;j<3;j++) Sedges.back().center[j] = face_center[j];
			}
			
			if (id2<0) //boundary faces
			{
				Bedges.push_back(edge(&Sgrains[id_curr],NULL,nbedges++,face_area));
				Bedges.back().active = id2;
				Sgrains[id_curr].bfaces.push_back(nbedges-1);
				
				std::string curr_string = entries[i+2*numfaces];
				curr_string = curr_string.substr(1,curr_string.size()-2);
				std::stringstream ss(curr_string);
				std::vector<int> face_vertex_ids;
				
				while( ss.good() )
				{
				    std::string substr;
				    std::getline(ss,substr,',');
				    face_vertex_ids.push_back(atoi(substr.c_str()));
				}
				
				for(int j=0;j<face_vertex_ids.size();j++)
				{
					//if(j<(face_vertex_ids.size()-1)) std::cout << "|";
					//else std::cout << " ";
					int local_vid = face_vertex_ids[j];
					curr_string = entries[local_vid+6+3*numfaces];
					curr_string = curr_string.substr(1,curr_string.size()-2);
					std::stringstream ss2(curr_string);
					
					Bedges.back().vertex_strs.push_back(entries[local_vid+6+3*numfaces]);
					
					while( ss2.good() )
					{
						std::string substr;
					    std::getline(ss2,substr,',');
					    
					    Bedges.back().vertices.push_back(::atof(substr.c_str()));
					    //result.push_back(::atof(substr.c_str()));
					}
					
					//std::cout << face_vertex_ids[j] << " [" << curr_string << "] ";
					//if(j<(face_vertex_ids.size()-1)) std::cout << "|";
					//else std::cout << " ";
						
				}
			}
			
			//if (id2>id_curr)
			//{
				//for(int j=0; j<Sedges.back().vertices.size(); j++)
				//{
					//std::cout << Sedges.back().vertices[j] << " ";
				//}
				//std::cout << std::endl;
			//}
			
			//if (id2==-1) leftF.push_back(id_curr);
			//if (id2==-2) rightF.push_back(id_curr);
		}
				
		
		//for(int i=0;i<entries.size();i++)
		//{
			//if(i!=0) outLfile << " ";
			//outLfile << entries[i];
		//}
		//outLfile << std::endl;
		
		//std::cout << numfaces << " " << numvertices << std::endl;
		
		//Sedges.push_back(edge(&Sgrains[i],&Sgrains[idn],nedges++,dd));
	}
	
	
	//std::cout << "solid MYSOLID";
	//std::cout << "Boundary edges = " << Bedges.size() << std::endl;
	
	
	//##// 2. Get face percolated structure via max-flow (after randomly activating voronoi faces, which represent graph edges)
	
	int perc_valid = 0,countS,num,num0;
	
	for(int i=0;i<Sedges.size();i++)
	{
		int m1 = Sedges[i].m->mol_id, n1 = Sedges[i].n->mol_id;
		if(m1!=n1) Sedges[i].active = 0;
		else Sedges[i].active = 1;
		//Sedges[i].active = 1;
		
		double upper = (COMS[(Sedges[i].m->mol_id)].z > COMS[(Sedges[i].n->mol_id)].z)? COMS[(Sedges[i].m->mol_id)].z : COMS[(Sedges[i].n->mol_id)].z;
		double lower = (COMS[(Sedges[i].m->mol_id)].z < COMS[(Sedges[i].n->mol_id)].z)? COMS[(Sedges[i].m->mol_id)].z : COMS[(Sedges[i].n->mol_id)].z;
		
		
		std::vector<int> cplane_id;
		for(int j=0;j<crackplanes.size();j++) if(lower<crackplanes[j] && upper>crackplanes[j]) cplane_id.push_back(j);
		
		//if(upper!=lower && Sedges[i].pbc_face!=1) 
		//if(upper!=lower && fabs(Sedges[i].center[2]-)<bounds.zhi )
		if(upper!=lower && ( (Sedges[i].pbc_face!=1)||(Sedges[i].pbc_face==1 && Sedges[i].pbc_zface==0 ) ) )
		{
			//printf("\n Sedge/face %d %4.12g %4.12g crack planes ",i,lower,upper);
			for(int j=0;j<cplane_id.size();j++) 
			{
				//printf(" %d %4.12g",cplane_id[j],crackplanes[cplane_id[j]]);
				if (Sedges[i].m->id<Sedges[i].n->id) 
				{
					crackplane_faces[cplane_id[j]].push_back(i);
					
					if (std::find(Sedges[i].plane_face.begin(), Sedges[i].plane_face.end(), cplane_id[j]) == Sedges[i].plane_face.end())
						Sedges[i].plane_face.push_back(cplane_id[j]);
				}
			}
		}
	}
	
	std::vector<double> areas(crackplanes.size(),0.0);
	std::vector<double> areasX(crackplanes.size(),0.0);
	std::vector<double> areasY(crackplanes.size(),0.0);
	std::vector<double> areasZ(crackplanes.size(),0.0);
	
	std::vector< std::set<point2> > projZ_points(crackplanes.size());
	std::map<point2, int> vertex_points;
	
	//printf("\n Areas");
	//for(int i=0;i<crackplanes.size();i++) printf(" %4.12g",areas[i]);
	
	for(int i=0;i<crackplane_faces.size();i++)
	{
		double min_z=bounds.zhi*2.0;
		
		for(int j=0;j<crackplane_faces[i].size();j++) 
		{
			int edge_id = crackplane_faces[i][j];
			
			double u[3] = {Sedges[edge_id].vertices[3]-Sedges[edge_id].vertices[0],Sedges[edge_id].vertices[4]-Sedges[edge_id].vertices[1],Sedges[edge_id].vertices[5]-Sedges[edge_id].vertices[2]}; 
			double v[3] = {Sedges[edge_id].vertices[6]-Sedges[edge_id].vertices[0],Sedges[edge_id].vertices[7]-Sedges[edge_id].vertices[1],Sedges[edge_id].vertices[8]-Sedges[edge_id].vertices[2]}; 
			double uxv[3] = {u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
			double nmag = sqrt(uxv[0]*uxv[0] + uxv[1]*uxv[1] + uxv[2]*uxv[2]);
			double n[3] = {uxv[0]/nmag,uxv[1]/nmag,uxv[2]/nmag};
			double cos = n[0]*direction[0] + n[1]*direction[1] + n[2]*direction[2];
			
			areasX[i] += Sedges[edge_id].len*fabs(n[0]);
			areasY[i] += Sedges[edge_id].len*fabs(n[1]);
			areasZ[i] += Sedges[edge_id].len*(n[2]);
			
			for(int k=0;k<Sedges[edge_id].vertices.size()/3;k++)
			{
				//if(Sedges[edge_id].pbc_face==1) 
				projZ_points[i].insert(point2(Sedges[edge_id].vertices[k*3],Sedges[edge_id].vertices[k*3+1],Sedges[edge_id].vertices[k*3+2]));
				if(min_z > Sedges[edge_id].vertices[k*3+2]) min_z = Sedges[edge_id].vertices[k*3+2];
			}
			
			areas[i] += Sedges[edge_id].len;
		}
		
		crackplanes_zmin[i] = min_z;
	}
	
	for(int i=0;i<crackplane_faces.size();i++)
	{
		std::stringstream buffer_temp;
		buffer_temp << "solids/Crackface" << i << "_points.dat";
		std::ofstream outfile_plane(buffer_temp.str().c_str());
		int vertex_count=0;
		
		for(std::set<point2>::const_iterator it = projZ_points[i].begin(); it != projZ_points[i].end(); ++it)
		{
			outfile_plane << it->x << " " << it->y << " " << it->z << "\n";
			
			if( fabs(it->z-crackplanes_zmin[i])<EPS2 )
			{
				if(vertex_points.find(*it) == vertex_points.end())
				{
					vertex_points.insert(std::pair<point2,int>(*it,vertex_count++));
				}
				
				//if(vertex_names.find(Sedges[curr_edge].vertex_strs[(order[i1])]) == vertex_names.end())
				//{
					//faces.push_back(vertex_count);
					//Svertices.push_back(point(Sedges[curr_edge].vertices[order[i1]*3],Sedges[curr_edge].vertices[order[i1]*3+1],Sedges[curr_edge].vertices[order[i1]*3+2],vertex_count,0));
					
					//Svertices[vertex_count].grains.insert(curr_grain);
					//vertex_names.insert(std::pair<std::string,int>(Sedges[curr_edge].vertex_strs[(order[i1])],vertex_count++));
				//}
				//else 
				//{
					//int vloc = vertex_names.at(   Sedges[curr_edge].vertex_strs[(order[i1])]   );
					//faces.push_back(vloc);
					//Svertices[vloc].grains.insert(curr_grain);
				//}
			}
			//std::cout << " " << *it;
			//int grain_id = *it;
			//center[0] += Sgrains[grain_id].x; center[1] += Sgrains[grain_id].y; center[2] += Sgrains[grain_id].z;
		}
		outfile_plane.close();
	}	
	
	double base_area = (bounds.xhi-bounds.xlo)*(bounds.yhi-bounds.ylo);
	
	printf("\n\nCrackface_id Total_area AreaX frac AreaY frac AreaZ* frac min_z");	
	for(int i=0;i<crackplanes.size();i++) printf("\n%d %4.12g %4.12g %4.12g %4.12g %4.12g %4.12g %4.12g %4.12g",i,areas[i],areasX[i],areasX[i]/base_area,areasY[i],areasY[i]/base_area,fabs(areasZ[i]),fabs(areasZ[i]/base_area),crackplanes_zmin[i]);
	//printf("\n Projected AreaX: ");	for(int i=0;i<crackplanes.size();i++) printf(" %4.12g",areasX[i]);
	//printf("\n Projected AreaY: ");	for(int i=0;i<crackplanes.size();i++) printf(" %4.12g",areasY[i]);
	//printf("\n Projected AreaZ: ");	for(int i=0;i<crackplanes.size();i++) printf(" %4.12g",fabs(areasZ[i]));
	//printf("\n Projected Area-fractions: ");	for(int i=0;i<crackplanes.size();i++) printf(" %4.12g",areasp[i]/areas[i]);
	
	using namespace boost;
	{
		
		//##// 3. Create new graph of percolated voronoi network, with nodes for voronoi cells and edges for voronoi faces
		
		typedef boost::property<boost::vertex_name_t, std::string> NameProperty;
		typedef boost::adjacency_list < boost::listS, boost::vecS, boost::undirectedS,
		    NameProperty > Graph;
		 
		typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;
		 
		typedef boost::property_map < Graph, boost::vertex_index_t >::type IndexMap;
		typedef boost::property_map < Graph, boost::vertex_name_t >::type NameMap;
			
		Graph G;
		std::vector <Vertex> vertices;
		std::map<int,int> vmap;
    
		for(int i=0;i<Sgrains.size();i++)
		{
			std::stringstream ss; ss << Sgrains[i].id;
			vertices.push_back(boost::add_vertex(ss.str(),G));
		}
		
		for(int i=0;i<Sedges.size();i++)
		{
			if(Sedges[i].active) 
			{
				add_edge(vertices[Sedges[i].m->id], vertices[Sedges[i].n->id], G);
			}//add_edge(Sedgest[i].m->id, Sedgest[i].n->id, G);
		}
	    
		std::vector<int> component(num_vertices(G));
		num = connected_components(G, &component[0]);
	    
		//std::cout << "Components " << num << std::endl;
    
		int count_bonds = 0;
		for(int i=0;i<Sedges.size();i++)
	    	//if(component[Sedges[i].m->id]==component[Sedges[i].n->id])
		if(Sedges[i].active) count_bonds++;
				
		std::cout << "\nTot_Sedges " << Sedges.size() << " bonds " << count_bonds << " frac " << count_bonds*1.0/Sedges.size() << " components " << num << std::endl;
		
		std::vector<complist> comp(num);
	    for(int i=0;i<Sgrains.size();i++)
	    {
				int compid = component[i];
				comp[compid].contents.push_back(i);
		}
		
		
		//int max=0,comp_max=0;
		//for(int i=0;i<comp.size();i++)
		//{
			//std::cout << std::endl << " Component " << i << " has " << comp[i].contents.size() << " grains :";
			//if( comp[i].contents.size() >= max)
			//{
				//max = comp[i].contents.size();
				//comp_max = i;
				////std::cout << "\n Hello\n";
			//}
			//for(int j=0;j<comp[i].contents.size();j++)
			//{
				//std::cout << " " << comp[i].contents[j];
			//}
		//}
		
		
		//std::cout << std::endl << "Largest Component " << comp_max << " which has " << max << " grains\n\n";
		
		
		//##// 4. Make Histogram
		
		std::map<int, int> comp_sizes;
		for(int i0=0;i0<(num);i0+=1)
		{
			int size0 = comp[i0].contents.size();
			if(comp_sizes.find(size0) == comp_sizes.end())
				comp_sizes.insert(std::pair<int,int>(size0,1));
			else
				comp_sizes.find(size0)->second++;
			
		}
		
		std::cout << "\nHistogram";
		for(std::map<int,int>::iterator it = comp_sizes.begin(); it != comp_sizes.end(); ++it)
		{
			std::cout << std::endl << it->first << " " << it->second;
		}
		std::cout << "\n";
		
		std::stringstream buffer;
		buffer << "solids/Assembly_all.stl";
		std::ofstream outfile(buffer.str().c_str());
		
		std::vector<std::ofstream> outfilec(crackplane_faces.size());
		std::vector<std::ofstream> outfilec2(crackplane_faces.size());
		std::vector<std::ofstream> outfilec3(crackplane_faces.size());
		
		for(int i=0;i<crackplane_faces.size();i++)
		{
			std::stringstream bufferc,bufferc2,bufferc3;
			bufferc << "solids/Assembly_allcracks_" << i << ".stl";
			outfilec[i].open(bufferc.str().c_str());
			outfilec[i] << "solid \"MYSOLID_planes\"";
			
			bufferc2 << "solids/Solid_below_crackplane_" << i << ".stl";
			outfilec2[i].open(bufferc2.str().c_str());
			outfilec2[i] << "solid \"MYSOLID\"";
			
			bufferc3 << "solids/Solid_above_crackplane_" << i << ".stl";
			outfilec3[i].open(bufferc3.str().c_str());
			outfilec3[i] << "solid \"MYSOLID\"";
			//std::ofstream out(bufferc.str().c_str());
			//out << "solid \"MYSOLID_planes\"";
			
			//outfilec.push_back(std::move(out));
		}
		//##// 5. Create STL file containing all components
		
		//for(int i0=0;i0<(num);i0+=1)
		for(int i=0;i<(num);i+=1)
		//for(int i0=comp_max;i0<=(comp_max+1);i0+=1)
		{
				int vertex_count = 0;
				
				std::vector<point> Svertices; //stores all the bounding vertices and the grains to which each vertex belongs to
				std::map<std::string, int> vertex_names;
				std::vector<int> faces;
				
				// SCANNING face vertices and storing vertex grain-affiliations for gap-displacement
				for(int j=0;j<comp[i].contents.size();j++)
				{
					int curr_grain = comp[i].contents[j];
					int grain2;
					//looping over internal faces
					for(int k=0;k<Sgrains[curr_grain].sfaces.size();k++)
					{
						int curr_edge = Sgrains[curr_grain].sfaces[k];
						
						if(Sedges[curr_edge].m->id == curr_grain)
							grain2 = Sedges[curr_edge].n->id;
						else if(Sedges[curr_edge].n->id == curr_grain)
							grain2 = Sedges[curr_edge].m->id;
						
						int cid2 = component[grain2];
						
						double dx0=Sgrains[grain2].x-Sgrains[curr_grain].x, dy0=Sgrains[grain2].y-Sgrains[curr_grain].y, dz0=Sgrains[grain2].z-Sgrains[curr_grain].z;
						minimum_image(dx0, dy0, dz0, bounds);
											
						//if(  (cid2!=i) || ( (cid2==i)&&(curr_grain<grain2) ) || ( (fabs(dx0)<EPS)&&(fabs(dy0)<EPS)&&(fabs(dz0)<EPS) )  ) // this ensures only faces across different components are treated AND internal faces are treated exactly only once for the smaller internal component
						//if(  (cid2!=i)
						if(  (cid2!=i) || (Sedges[curr_edge].pbc_face==1) ) // this ensures only faces across different components are treated AND pbc faces too
						{
							double dx=Sgrains[grain2].x-Sgrains[curr_grain].x, dy=Sgrains[grain2].y-Sgrains[curr_grain].y, dz=Sgrains[grain2].z-Sgrains[curr_grain].z;
							minimum_image(dx, dy, dz, bounds);
							//double dir[3] = {dx,dy,dz};
							double dir[3] = {Sedges[curr_edge].normal[0],Sedges[curr_edge].normal[1],Sedges[curr_edge].normal[2]};
							
							//std::cout << "\n " << curr_grain << " " << grain2 << " " << " displacement " << dx << " " << dy << " " << dz;
							
							//double dir[3] = {Sgrains[grain2].x-Sgrains[curr_grain].x,Sgrains[grain2].y-Sgrains[curr_grain].y,Sgrains[grain2].z-Sgrains[curr_grain].z};
							double dirmag = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
							double dirn[3] = {dir[0]/dirmag,dir[1]/dirmag,dir[2]/dirmag};
							
							double u[3] = {Sedges[curr_edge].vertices[3]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[4]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[5]-Sedges[curr_edge].vertices[2]}; 
							double v[3] = {Sedges[curr_edge].vertices[6]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[7]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[8]-Sedges[curr_edge].vertices[2]}; 
							double uxv[3] = {u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
							double nmag = sqrt(uxv[0]*uxv[0] + uxv[1]*uxv[1] + uxv[2]*uxv[2]);
							double n[3] = {uxv[0]/nmag,uxv[1]/nmag,uxv[2]/nmag};
							
							double dot = dirn[0]*n[0]+dirn[1]*n[1]+dirn[2]*n[2];
							
							if(fabs(dirmag)<EPS) {dirn[0] = n[0];dirn[1] = n[1];dirn[2] = n[2];dot=1;}
							//outfile << std::endl << "Dot " << dot;
							
							//int order[3] = {0,1,2};
							//if(dot<0) {order[1]=2; order[2]=1;}
							
							//if(fabs(dot)<EPS) std::cout << "\n " << curr_grain << " " << grain2 << " " << " displacement " << dx << " " << dy << " " << dz;
							
							int num_edge_verts = Sedges[curr_edge].vertices.size()/3;
							
							//std::cout << "\n Num_edge_vertices " << num_edge_verts;
							
							//if (Sedges[curr_edge].vertices.size()==9)
							for(int i0=0;i0<(num_edge_verts-2);i0++)
							{
								int a0,b0;
								if(dot>=0) {a0 = i0+1; b0=i0+2;}
								else {a0 = i0+2; b0=i0+1;}
								
								int order[3] = {0,a0,b0};
								
								for(int i1=0; i1<3; i1++)
								{
									if(vertex_names.find(Sedges[curr_edge].vertex_strs[(order[i1])]) == vertex_names.end())
									{
										faces.push_back(vertex_count);
										Svertices.push_back(point(Sedges[curr_edge].vertices[order[i1]*3],Sedges[curr_edge].vertices[order[i1]*3+1],Sedges[curr_edge].vertices[order[i1]*3+2],vertex_count,0));
										
										Svertices[vertex_count].grains.insert(curr_grain);
										vertex_names.insert(std::pair<std::string,int>(Sedges[curr_edge].vertex_strs[(order[i1])],vertex_count++));
									}
									else 
									{
										int vloc = vertex_names.at(   Sedges[curr_edge].vertex_strs[(order[i1])]   );
										faces.push_back(vloc);
										Svertices[vloc].grains.insert(curr_grain);
									}
								}
							}
						}
					}
					
					//looping over boundary faces, internal component should have no boundary faces
					for(int k=0;k<Sgrains[curr_grain].bfaces.size();k++)
					{
						int curr_edge = Sgrains[curr_grain].bfaces[k];
						int facen = -Bedges[curr_edge].active;
						
						double dirn[3] = {0,0,0};
						dirn[(facen-1)/2] = 1-2*(facen%2);
						double grn[3] = {0,0,0}; 
						double toln[3] = {0,0,0};
						
						
						double u[3] = {Bedges[curr_edge].vertices[3]-Bedges[curr_edge].vertices[0],Bedges[curr_edge].vertices[4]-Bedges[curr_edge].vertices[1],Bedges[curr_edge].vertices[5]-Bedges[curr_edge].vertices[2]}; 
						double v[3] = {Bedges[curr_edge].vertices[6]-Bedges[curr_edge].vertices[0],Bedges[curr_edge].vertices[7]-Bedges[curr_edge].vertices[1],Bedges[curr_edge].vertices[8]-Bedges[curr_edge].vertices[2]}; 
						double uxv[3] = {u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
						double nmag = sqrt(uxv[0]*uxv[0] + uxv[1]*uxv[1] + uxv[2]*uxv[2]);
						double n[3] = {uxv[0]/nmag,uxv[1]/nmag,uxv[2]/nmag};
						
						double dot = dirn[0]*n[0]+dirn[1]*n[1]+dirn[2]*n[2],nmag2;
						
						//int order[3] = {0,1,2};
						//if(dot<0) {order[1]=2; order[2]=1;}
						
						int num_edge_verts = Bedges[curr_edge].vertices.size()/3;
						
						for(int i0=0;i0<(num_edge_verts-2);i0++)
						{
							int a0,b0;
							if(dot>=0) {a0 = i0+1; b0=i0+2;}
							else {a0 = i0+2; b0=i0+1;}
							
							int order[3] = {0,a0,b0};
							
							for(int i1=0; i1<3; i1++)
							{
								if(vertex_names.find(Bedges[curr_edge].vertex_strs[(order[i1])]) == vertex_names.end())
								{
									faces.push_back(vertex_count);
									Svertices.push_back(point(Bedges[curr_edge].vertices[order[i1]*3],Bedges[curr_edge].vertices[order[i1]*3+1],Bedges[curr_edge].vertices[order[i1]*3+2],vertex_count,0));
									
									Svertices[vertex_count].grains.insert(curr_grain);
									vertex_names.insert(std::pair<std::string,int>(Bedges[curr_edge].vertex_strs[(order[i1])],vertex_count++));
								}
								else 
								{
									int vloc = vertex_names.at(   Bedges[curr_edge].vertex_strs[(order[i1])]   );
									faces.push_back(vloc);
									Svertices[vloc].grains.insert(curr_grain);
								}
							}
						}
					}
				}
				
				//if(i==0) outfile << "solid \"MYSOLID_" << i << "." << component_i[comp[i0].contents[0]] << "\"";
				//else outfile << "solid \"MYSOLID_" << i << "\"";
				//outfile1 << "solid \"MYSOLID_" << i << "\"";
				
				std::stringstream buffer_temp;
				buffer_temp << "solids/MYSOLID" << i << ".stl";
				std::ofstream outfile_solid(buffer_temp.str().c_str());
				
				outfile << "solid \"MYSOLID_" << i << "\"";
				outfile_solid << "solid \"MYSOLID_" << i << "\"";
				
				/////////////////  printing STL faces //////////////
				for(int j=0;j<comp[i].contents.size();j++)
				{
					int curr_grain = comp[i].contents[j];
					int grain2;
					
					//std::cout << "\n comp " << i << " has faces " << Sgrains[curr_grain].sfaces.size() << " " << Sgrains[curr_grain].bfaces.size();
					
					//looping over internal faces
					for(int k=0;k<Sgrains[curr_grain].sfaces.size();k++)
					{
						int curr_edge = Sgrains[curr_grain].sfaces[k];
						
						if(Sedges[curr_edge].m->id == curr_grain)
							grain2 = Sedges[curr_edge].n->id;
						else if(Sedges[curr_edge].n->id == curr_grain)
							grain2 = Sedges[curr_edge].m->id;
						
						int cid2 = component[grain2];
						
						double dx0=Sgrains[grain2].x-Sgrains[curr_grain].x, dy0=Sgrains[grain2].y-Sgrains[curr_grain].y, dz0=Sgrains[grain2].z-Sgrains[curr_grain].z;
						minimum_image(dx0, dy0, dz0, bounds);
						
						//if(cid2 != i) * removed because in the inclusion all internal faces must be present so that crack faces in the .cracks can be matched with mesh faces that have been previously defined
						//if(  (cid2!=i) || ( (cid2==i)&&(curr_grain<grain2) )  ) // this ensures only faces across different components are treated AND internal faces are treated exactly only once
						//if(  (cid2!=i) || ( (cid2==i)&&(curr_grain<grain2) ) || ( (fabs(dx0)<EPS)&&(fabs(dy0)<EPS)&&(fabs(dz0)<EPS) )  )
						if(  (cid2!=i) || (Sedges[curr_edge].pbc_face==1) ) // this ensures only faces across different components are treated AND pbc faces too
						{
							double dx=Sgrains[grain2].x-Sgrains[curr_grain].x, dy=Sgrains[grain2].y-Sgrains[curr_grain].y, dz=Sgrains[grain2].z-Sgrains[curr_grain].z;
							minimum_image(dx, dy, dz, bounds);
							//double dir[3] = {dx,dy,dz};
							double dir[3] = {Sedges[curr_edge].normal[0],Sedges[curr_edge].normal[1],Sedges[curr_edge].normal[2]};
							
							//double dir[3] = {Sgrains[grain2].x-Sgrains[curr_grain].x,Sgrains[grain2].y-Sgrains[curr_grain].y,Sgrains[grain2].z-Sgrains[curr_grain].z};
							double dirmag = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
							double dirn[3] = {dir[0]/dirmag,dir[1]/dirmag,dir[2]/dirmag};
							
							double u[3] = {Sedges[curr_edge].vertices[3]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[4]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[5]-Sedges[curr_edge].vertices[2]}; 
							double v[3] = {Sedges[curr_edge].vertices[6]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[7]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[8]-Sedges[curr_edge].vertices[2]}; 
							double uxv[3] = {u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
							double nmag = sqrt(uxv[0]*uxv[0] + uxv[1]*uxv[1] + uxv[2]*uxv[2]);
							double n[3] = {uxv[0]/nmag,uxv[1]/nmag,uxv[2]/nmag};
							
							double dot = dirn[0]*n[0]+dirn[1]*n[1]+dirn[2]*n[2]; // dot=1 or -1
							//outfile << std::endl << "Dot " << dot;
							
							//std::cout << "\n dirmag " << dirmag << " " << dot;
							if(fabs(dirmag)<EPS) std::cout << "\n dirmag " << dirmag << " " << dot;
							
							int sign = (dot > 0) - (dot < 0);
							dirn[0] = sign*n[0];dirn[1] = sign*n[1];dirn[2] = sign*n[2];
						
							Sedges[curr_edge].normal2[0] = dirn[0];
							//if(fabs(dot)<EPS) std::cout << "\n " << curr_grain << " " << grain2 << " " << " displacement " << dx << " " << dy << " " << dz;
							
							//int order[3] = {0,1,2};
							//if(dot<0) {order[1]=2; order[2]=1;}
							
							int num_edge_verts = Sedges[curr_edge].vertices.size()/3;
														
							//if (Sedges[curr_edge].vertices.size()==12) // EXTRA triangular FACE if rhombus
							for(int i0=0;i0<(num_edge_verts-2);i0++)
							{
								int a0,b0;
								if(dot>=0) {a0 = i0+1; b0=i0+2;}
								else {a0 = i0+2; b0=i0+1;}
								
								int order[3] = {0,a0,b0};
								
								//order[1] += 1; order[2] += 1;
								
								//if(cid2!=i) 
								outfile << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
								outfile_solid << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
								
								//if((Sedges[curr_edge].plane_face!=-1)) 
								for(int ip=0;ip<Sedges[curr_edge].plane_face.size();ip++)
									outfilec[ Sedges[curr_edge].plane_face[ip] ] << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
								
								for(int ip=0;ip<crackplane_faces.size();ip++)
								{
									if (COMS[(Sedges[curr_edge].m->mol_id)].z<crackplanes[ip] && (COMS[(Sedges[curr_edge].n->mol_id)].z>crackplanes[ip] || Sedges[curr_edge].pbc_face==1 ))
										outfilec2[ip] << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
									
									//if (COMS[(Sedges[curr_edge].m->mol_id)].z>crackplanes[ip])
									int molid1 = (ip/2+1)*nx*ny*2+(z_switch+2*(nx*ny/2)), molid2 = ((ip+1)/2)*nx*ny*2+(1-z_switch+2*(nx*ny/2));
									if( (Sedges[curr_edge].m->mol_id==molid1 && Sedges[curr_edge].n->mol_id!=molid2 )
											|| (Sedges[curr_edge].m->mol_id==molid2 && Sedges[curr_edge].n->mol_id!=molid1) )
									outfilec3[ip] << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
								}
								//COMS[(Sedges[i].m->mol_id)].z
								
								//*if((Sedges[curr_edge].plane_face==1)) outfile2 << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
								//outfile1 << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
																
								for(int i1=0; i1<3; i1++)
								{
									int vloc = vertex_names.at(   Sedges[curr_edge].vertex_strs[(order[i1])]   );
									//std::set<int>::iterator it;
									double center[] = {0,0,0};
									
									//std::cout << "\n" << Svertices[vloc].id << "Vertex " << Sedges[curr_edge].vertex_strs[(order[i1])] << " " << Svertices[vloc].id << " " << Svertices[vloc].x << " " << Svertices[vloc].y << " " << Svertices[vloc].z << " (" << Svertices[vloc].grains.size() << ") ";
									
									for(std::set<int>::const_iterator it = Svertices[vloc].grains.begin(); it != Svertices[vloc].grains.end(); ++it)
									{
										//std::cout << " " << *it;
										int grain_id = *it;
										center[0] += Sgrains[grain_id].x; center[1] += Sgrains[grain_id].y; center[2] += Sgrains[grain_id].z;
									}
									for(int i2=0;i2<3;i2++) center[i2] = center[i2]/Svertices[vloc].grains.size();
									
									double grn[] = { Svertices[vloc].x - center[0], Svertices[vloc].y - center[1], Svertices[vloc].z - center[2] };
									double nmag2 = sqrt(grn[0]*grn[0]+grn[1]*grn[1]+grn[2]*grn[2]); 
									double toln[] = { grn[0]/nmag2, grn[1]/nmag2, grn[2]/nmag2 };
									if (nmag2<EPS) {toln[0]=0.0; toln[1]=0.0; toln[2]=0.0;}
									
									//if(cid2!=i) 
									outfile << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
									outfile_solid << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
									
									//if((Sedges[curr_edge].plane_face!=-1)) 
									for(int ip=0;ip<Sedges[curr_edge].plane_face.size();ip++)
										outfilec[ Sedges[curr_edge].plane_face[ip] ] << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
									
									for(int ip=0;ip<crackplane_faces.size();ip++)
									{
										if (COMS[(Sedges[curr_edge].m->mol_id)].z<crackplanes[ip] && (COMS[(Sedges[curr_edge].n->mol_id)].z>crackplanes[ip] || Sedges[curr_edge].pbc_face==1 ))
											outfilec2[ip] << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
										
										//if (COMS[(Sedges[curr_edge].m->mol_id)].z>crackplanes[ip])
										//int molid1 = (ip/2+1)*nx*ny*2+8, molid2 = ((ip+1)/2)*nx*ny*2+9;
										int molid1 = (ip/2+1)*nx*ny*2+(z_switch+2*(nx*ny/2)), molid2 = ((ip+1)/2)*nx*ny*2+(1-z_switch+2*(nx*ny/2));
										//if( Sedges[curr_edge].m->mol_id==molid1 || Sedges[curr_edge].m->mol_id==molid2 )
										if( (Sedges[curr_edge].m->mol_id==molid1 && Sedges[curr_edge].n->mol_id!=molid2 )
											|| (Sedges[curr_edge].m->mol_id==molid2 && Sedges[curr_edge].n->mol_id!=molid1) )
											outfilec3[ip] << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
									}
									//*if((Sedges[curr_edge].plane_face==1)) outfile2 << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
								
								
									//if((Sedges[curr_edge].plane_face==1)) outfile2 << "(" << Svertices[vloc].x-TOL*toln[0] << "," << Svertices[vloc].y-TOL*toln[1] << "," << Svertices[vloc].z-TOL*toln[2] << ") ";
									//if((Sedges[curr_edge].active==0)&&(cid2==i)) outfile2 << "(" << Svertices[vloc].x-TOL*toln[0] << "," << Svertices[vloc].y-TOL*toln[1] << "," << Svertices[vloc].z-TOL*toln[2] << ") ";
									//outfile1 << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2] << " " << Svertices[vloc].id << "Vertex " << Sedges[curr_edge].vertex_strs[(order[i1])] << " " << Svertices[vloc].id << " " << Svertices[vloc].x << " " << Svertices[vloc].y << " " << Svertices[vloc].z << " (" << Svertices[vloc].grains.size() << ") ";
									//std::cout << std::endl << Svertices[vloc].id << "      Fvertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
								}
								
								//if(cid2!=i) 
								outfile << std::endl << "    endloop" << std::endl << "  endfacet";
								outfile_solid << std::endl << "    endloop" << std::endl << "  endfacet";
								
								for(int ip=0;ip<Sedges[curr_edge].plane_face.size();ip++)
									outfilec[Sedges[curr_edge].plane_face[ip] ] << std::endl << "    endloop" << std::endl << "  endfacet";
								//*if((Sedges[curr_edge].plane_face==1)) outfile2 << std::endl << "    endloop" << std::endl << "  endfacet";
								
								for(int ip=0;ip<crackplane_faces.size();ip++)
								{
									if (COMS[(Sedges[curr_edge].m->mol_id)].z<crackplanes[ip] && (COMS[(Sedges[curr_edge].n->mol_id)].z>crackplanes[ip] || Sedges[curr_edge].pbc_face==1 ))
										outfilec2[ip] << std::endl << "    endloop" << std::endl << "  endfacet";
										
									//if (COMS[(Sedges[curr_edge].m->mol_id)].z>crackplanes[ip])
									//int molid1 = (ip/2+1)*nx*ny*2+8, molid2 = ((ip+1)/2)*nx*ny*2+9;
									int molid1 = (ip/2+1)*nx*ny*2+(z_switch+2*(nx*ny/2)), molid2 = ((ip+1)/2)*nx*ny*2+(1-z_switch+2*(nx*ny/2));
									//if( Sedges[curr_edge].m->mol_id==molid1 || Sedges[curr_edge].m->mol_id==molid2 )
									if( (Sedges[curr_edge].m->mol_id==molid1 && Sedges[curr_edge].n->mol_id!=molid2 )
											|| (Sedges[curr_edge].m->mol_id==molid2 && Sedges[curr_edge].n->mol_id!=molid1) )
									outfilec3[ip] << std::endl << "    endloop" << std::endl << "  endfacet";	
								}
								//if((Sedges[curr_edge].plane_face==1)) outfile2 << std::endl;
								//outfile1 << std::endl << "    endloop" << std::endl << "  endfacet";
								
							}
						}
					}
					
					//looping over boundary faces
					for(int k=0;k<Sgrains[curr_grain].bfaces.size();k++)
					{
						int curr_edge = Sgrains[curr_grain].bfaces[k];
						int facen = -Bedges[curr_edge].active;
						
						double dirn[3] = {0,0,0};
						dirn[(facen-1)/2] = 1-2*(facen%2);
						double grn[3] = {0,0,0}; 
						double toln[3] = {0,0,0};
						
						
						double u[3] = {Bedges[curr_edge].vertices[3]-Bedges[curr_edge].vertices[0],Bedges[curr_edge].vertices[4]-Bedges[curr_edge].vertices[1],Bedges[curr_edge].vertices[5]-Bedges[curr_edge].vertices[2]}; 
						double v[3] = {Bedges[curr_edge].vertices[6]-Bedges[curr_edge].vertices[0],Bedges[curr_edge].vertices[7]-Bedges[curr_edge].vertices[1],Bedges[curr_edge].vertices[8]-Bedges[curr_edge].vertices[2]}; 
						double uxv[3] = {u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
						double nmag = sqrt(uxv[0]*uxv[0] + uxv[1]*uxv[1] + uxv[2]*uxv[2]);
						double n[3] = {uxv[0]/nmag,uxv[1]/nmag,uxv[2]/nmag};
						
						double dot = dirn[0]*n[0]+dirn[1]*n[1]+dirn[2]*n[2],nmag2;
						
						//int order[3] = {0,1,2};
						//if(dot<0) {order[1]=2; order[2]=1;}
						
						int num_edge_verts = Bedges[curr_edge].vertices.size()/3;
														
						//if (Bedges[curr_edge].vertices.size()==12) // EXTRA triangular FACE if rhombus
						for(int i0=0;i0<(num_edge_verts-2);i0++)
						{
							int a0,b0;
							if(dot>=0) {a0 = i0+1; b0=i0+2;}
							else {a0 = i0+2; b0=i0+1;}
							
							int order[3] = {0,a0,b0};
							//order[1] += 1; order[2] += 1;
							
							outfile << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
							outfile_solid << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
							//outfile1 << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
							
							for(int i1=0; i1<3; i1++)
							{
								int vloc = vertex_names.at(   Bedges[curr_edge].vertex_strs[(order[i1])]   );
								std::set<int>::iterator it;
								double center[] = {0,0,0};
								
								for( it = Svertices[vloc].grains.begin(); it != Svertices[vloc].grains.end(); ++it)
								{
									int grain_id = *it;
									center[0] += Sgrains[grain_id].x; center[1] += Sgrains[grain_id].y; center[2] += Sgrains[grain_id].z;
								}
								
								for(int i2=0;i2<3;i2++) center[i2] = center[i2]/Svertices[vloc].grains.size();
								
								double grn[] = { Svertices[vloc].x - center[0], Svertices[vloc].y - center[1], Svertices[vloc].z - center[2] };
								double nmag2 = sqrt(grn[0]*grn[0]+grn[1]*grn[1]+grn[2]*grn[2]); 
								double toln[] = { grn[0]/nmag2, grn[1]/nmag2, grn[2]/nmag2 };
								if (nmag2<EPS) {toln[0]=0.0; toln[1]=0.0; toln[2]=0.0;}
								
								outfile << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
								outfile_solid << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
								
							}
							
							outfile << std::endl << "    endloop" << std::endl << "  endfacet";
							outfile_solid << std::endl << "    endloop" << std::endl << "  endfacet";
							//outfile1 << std::endl << "    endloop" << std::endl << "  endfacet";
						}
						
						
						
						////
					
					}
				
							
				}
				
				outfile << std::endl << "endsolid" << std::endl;
				outfile_solid << std::endl << "endsolid" << std::endl;
			
				outfile_solid.close();
				///////////////////////
				
				//for(std::map<std::string, int>::const_iterator it = vertex_names.begin(); it != vertex_names.end(); ++it)
				//{
					//std::cout << it->first << " " << it->second << " " << Svertices[it->second].x << " " << Svertices[it->second].y << " " << Svertices[it->second].z << "\n";
				//}
				
				//outfile << std::endl << "endsolid" << std::endl;
				//outfile1 << std::endl << "endsolid" << std::endl;
			
		}
		
		//*outfile2 << std::endl << "endsolid" << std::endl;
		
		outfile.close();
		
		for(int i=0;i<crackplane_faces.size();i++)
		{
			
			outfilec[i] << std::endl << "endsolid" << std::endl;
			outfilec[i].close();
			
			outfilec2[i] << std::endl << "endsolid" << std::endl;
			outfilec2[i].close();
			
			outfilec3[i] << std::endl << "endsolid" << std::endl;
			outfilec3[i].close();
			
		}
		//*outfile2.close();
			
		
		//for(std::map<std::string, int>::const_iterator it = vertex_names.begin(); it != vertex_names.end(); ++it)
		//{
			//std::cout << it->first << " " << it->second << " " << Svertices[it->second].x << " " << Svertices[it->second].y << " " << Svertices[it->second].z << "\n";
		//}
		
		//for(int i=0;i<Svertices.size();i++)
		//{
				//outfile << Svertices[i].x << " " << Svertices[i].y << " " << Svertices[i].z << "\n";
		//}
		
		//for(int i=0;i<faces.size();i++)
		//{
				//outfile1 << faces[i];
				//if(i%3==2) outfile1 << std::endl;
				//else outfile1 << " ";
		//}
		
		//outfile1.close();
		//outfile.close();
	    
	}
	
	std::cout << std::endl;
	
	//infile.close();
	////outfile.close();
	
	int stop_s=clock();
	std::cout << "Total time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << std::endl;
	
	return 0;
	
}

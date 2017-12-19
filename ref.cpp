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

#define TOL2 1e-6
#define MIN 2.0
#define MAX 3.0
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
  int cell_id;
  std::vector<int> sfaces,bfaces,edges;
  std::set<int> grains;
  point(double a, double b, double c, int d, int e) : x(a), y(b), z(c), id(d), cell_id(e) {}
};

struct edge {
  //double x,y,z;
  point *m,*n;
  int id,p,q,active;
  double len;
  std::vector<double> vertices;
  std::vector<std::string> vertex_strs;
  edge(point *a, point *b, int d, double l) : m(a), n(b), id(d), len(l) {}
};

struct complist {
	std::vector<int> contents;
};

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
		printf("\n USAGE: ./create_assembly_stl_lmp <fraction> <size> <percolated yes/no(1/0)>"); 
		exit(EXIT_FAILURE);
	}
	
	std::string line,thrash,thrash1;
	double frac0 = atof(argv[1]);
	double sizeN = lata*atof(argv[2]);
	int perc = atoi(argv[3]);
		
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
		
	std::ofstream outLfile("ABA.temp");
	std::ofstream outLfile2("data.lmp");
	std::ofstream outLfileG("groups.dat");
	std::stringstream sstream;
	
	std::vector<point> Sgrains;
	std::vector<edge> Sedges,Bedges;
	std::vector<int> leftF,rightF;
	
	for(int i=0;i<ngrains;i++) Sgrains.push_back(point(0.00001,0.00001,0.00001,i+1,0));
	
	//for(int i=0;i<ngrains;i++) std::cout << std::endl << Sgrains[i].id << " " << Sgrains[i].x << " " << Sgrains[i].y << " " << Sgrains[i].z;
			
	std::string input;
	
	infile0.close(); std::ifstream infile("voro_cell.dat");
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
		
		int id_curr = atoi(entries[0].c_str());
		Sgrains[id_curr].id = id_curr;
		Sgrains[id_curr].x = ::atof(entries[1].c_str()); Sgrains[id_curr].y = ::atof(entries[2].c_str()); Sgrains[id_curr].z = ::atof(entries[3].c_str());
		//std::cout << std::endl << id_curr << " " << Sgrains[id_curr].id << " " << Sgrains[id_curr].x << " " << Sgrains[id_curr].y << " " << Sgrains[id_curr].z;
		
		
		int numfaces = atoi(entries[4].c_str());
		for(int i=5;i<(5+numfaces);i++)
		{
			int id2 = atoi(entries[i].c_str());
			double face_area = ::atof(entries[i+numfaces].c_str());
			
			if (id2>id_curr) 
			{
				Sedges.push_back(edge(&Sgrains[id_curr],&Sgrains[id2],nedges++,face_area));
				Sgrains[id_curr].sfaces.push_back(nedges-1);
				Sgrains[id2].sfaces.push_back(nedges-1);
				
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
					
					while( ss2.good() )
					{
						std::string substr;
					    std::getline(ss2,substr,',');
					    
					    Sedges.back().vertices.push_back(::atof(substr.c_str()));
					    //result.push_back(::atof(substr.c_str()));
					}
					
					//std::cout << face_vertex_ids[j] << " [" << curr_string << "] ";
					//if(j<(face_vertex_ids.size()-1)) std::cout << "|";
					//else std::cout << " ";
						
				}
			}
			
			if (id2<0) 
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
					
			/*
			if (id2>id_curr)
			{
				for(int j=0; j<Sedges.back().vertices.size(); j++)
				{
					std::cout << Sedges.back().vertices[j] << " ";
				}
				std::cout << std::endl;
			}*/
			
			if (id2==-1) leftF.push_back(id_curr);
			if (id2==-2) rightF.push_back(id_curr);
		}
				
		/*
		for(int i=0;i<entries.size();i++)
		{
			if(i!=0) outLfile << " ";
			outLfile << entries[i];
		}
		outLfile << std::endl;
		*/
		//std::cout << numfaces << " " << numvertices << std::endl;
		
		//Sedges.push_back(edge(&Sgrains[i],&Sgrains[idn],nedges++,dd));
	}
	
	//std::cout << "solid MYSOLID";
	//std::cout << "Boundary edges = " << Bedges.size() << std::endl;
	
	
	//##// 2. Get face percolated structure via max-flow (after randomly activating voronoi faces, which represent graph edges)
	
	int perc_valid = 0,countS,num,num0;
		
	do
	{	
		using namespace boost;
		{
			//typedef boost::property<boost::vertex_name_t, std::string> NameProperty;
			typedef boost::adjacency_list < vecS, vecS, directedS,
				  property < vertex_name_t, std::string,
				    property < vertex_index_t, int,
				      property < vertex_color_t, boost::default_color_type,
				        property < vertex_distance_t, double_t,
				          property < vertex_predecessor_t, Traits::edge_descriptor > > > > >,
				 
				  property < edge_capacity_t, EdgeWeightType,
				    property < edge_residual_capacity_t, EdgeWeightType,
				      property < edge_reverse_t, Traits::edge_descriptor > > > > Graph;
			 
			typedef Traits::vertex_descriptor Vertex;
						
			typedef boost::property_map < Graph, boost::vertex_color_t >::type ColorMap;					 
			//typedef boost::property_map < Graph, boost::vertex_index_t >::type IndexMap;
			//typedef boost::property_map < Graph, boost::vertex_name_t >::type NameMap;
			
		    Graph G;
		    property_map < Graph, edge_reverse_t >::type rev = get(edge_reverse, G);
		    std::set<int> intersect;
		    std::vector<Vertex> vertices;
		    
		    int leftvertex=-1, rightvertex=-1;
			
			for(int i=0;i<Sgrains.size();i++)
		    {
				//std::stringstream ss; ss << Sgrainst[i].id;
				vertices.push_back(boost::add_vertex(G));
			}
			
			vertices.push_back(boost::add_vertex(G)); leftvertex = Sgrains.size(); //source
			vertices.push_back(boost::add_vertex(G)); rightvertex = Sgrains.size()+1; //sink
					    
			countS=0;
			for(int i=0;i<Sedges.size();i++)
			{
				double area;
				
				int m1 = Sedges[i].m->id, n1 = Sedges[i].n->id;
				Sgrains[m1].edges.push_back(i);	Sgrains[n1].edges.push_back(i);
								
				double d[3] = {Sgrains[n1].x - Sgrains[m1].x, Sgrains[n1].y - Sgrains[m1].y, Sgrains[n1].z - Sgrains[m1].z};
				double cos = fabs(Sgrains[n1].x - Sgrains[m1].x)/sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
										
				double rand01 = rn();
				//std::cout <<"\nRandom " << rand01;
				
				if(rand01>frac0) 
				{
					area = 100000.0;
					Sedges[i].active = 1;
				}
				else 
				{
					area = Sedges[i].len*cos;
					Sedges[i].active = 0;
				}
				/*
				if( fabs(Sgrainst[n1].x-0) <= TOL && fabs(Sgrainst[m1].x-0) <= TOL ) area = area*0.5;
				if( fabs(Sgrainst[n1].x-1.0*(NX-1)) <= TOL && fabs(Sgrainst[m1].x-1.0*(NX-1)) <= TOL ) area = area*0.5;
				
				if( fabs(Sgrainst[n1].y-0) <= TOL && fabs(Sgrainst[m1].y-0) <= TOL ) area = area*0.5;
				if( fabs(Sgrainst[n1].y-1.0*(NY-1)) <= TOL && fabs(Sgrainst[m1].y-1.0*(NY-1)) <= TOL ) area = area*0.5;
										
				if( fabs(Sgrainst[n1].z-0) <= TOL && fabs(Sgrainst[m1].z-0) <= TOL ) area = area*0.5;
				if( fabs(Sgrainst[n1].z-1.0*(NZ-1)) <= TOL && fabs(Sgrainst[m1].z-1.0*(NZ-1)) <= TOL ) area = area*0.5;
				*/					
				
				AddEdge(vertices[m1], vertices[n1], rev, area, G);
			
				//add_edge(vertices[m1], vertices[n1], G);
				countS++;
				//add_edge(Sedgest[i].m->id, Sedgest[i].n->id, G);
				
			}
			
			//Add additional joined edges between source and leftface && sink and rightface
			for(int i=0;i<leftF.size();i++)
			{
				double area = 100000.0;
				int idn = leftF[i];
				AddEdge(vertices[leftvertex], vertices[idn], rev, area, G);
			}
			
			for(int i=0;i<rightF.size();i++)
			{
				double area = 100000.0;
				int idn = rightF[i];
				AddEdge(vertices[rightvertex], vertices[idn], rev, area, G);
			}
			
			EdgeWeightType flow = boykov_kolmogorov_max_flow(G, vertices[leftvertex], vertices[rightvertex]);
		    ColorMap color = get(vertex_color, G);
		    
		    for(int i=0;i<Sgrains.size();i++)
		    {
				Sgrains[i].color = color[i];
				//std::cout << "\n color:" << i << "= " << atoi(color[i]);
			}
			
			for(int i=0;i<Sedges.size();i++)
			{
				int m1 = Sedges[i].m->id, n1 = Sedges[i].n->id;
				if(Sgrains[m1].color==Sgrains[n1].color) Sedges[i].active=1;
				else Sedges[i].active=0;
			}
		    
		    //std::vector<int> component(num_vertices(G));
			//int num = connected_components(G, &component[0]);
			
			//std::cout << "Max flow is: " << flow << std::endl;
			
			if(flow<10.0) perc_valid++ ;
		}
		
	}while((perc==1)&&(perc_valid<1));
	
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
				
		std::cout << "\nTot_Sedges " << Sedges.size() << " bonds " << count_bonds << " frac " << count_bonds*1.0/Sedges.size() << " frac0 " << frac0 << std::endl;
		    
		///////////////////////// WRITING LAMMPS FILES (data.lmp & groups.dat) ////////////////////////////
	
    outLfile << "\nAtoms\n\n";
    
    std::vector<int> F[3],L[3];
    std::set<int> rigidF[3],rigidL[3];
    std::vector<complist> comp(num);
    for(int i=0;i<Sgrains.size();i++)
    {
			int compid = component[i];
			comp[compid].contents.push_back(i);
			//outLfile << Sgrains[i].id+1 << " " << component[i]+1 << " 1 " << Sgrains[i].sfaces.size() << " " << (Sgrains[i].x-MIN)*sizeN << " " << (Sgrains[i].y-MIN)*sizeN << " " << (Sgrains[i].z-MIN)*sizeN << std::endl;
			outLfile << Sgrains[i].id+1 << " " << component[i]+1 << " 1 " << 0.0 << " " << std::round((Sgrains[i].x-MIN)*sizeN) << " " << std::round((Sgrains[i].y-MIN)*sizeN) << " " << std::round((Sgrains[i].z-MIN)*sizeN) << std::endl;
			std::cout << std::endl << i+1 << "-";
			//for(int j=0;j<Sgrains[i].sfaces.size();j++) std::cout << " " << Sgrains[i].sfaces[j];
			//for(int j=0;j<Sgrains[i].sfaces.size();j++)
			//{
				//int iter = Sgrains[i].sfaces[j];
				//int other_id = (Sedges[iter].m->id == i) ? Sedges[iter].n->id : Sedges[iter].m->id;
				//std::cout << " " << other_id+1 << " (" << Sedges[iter].m->id+1 << "," << Sedges[iter].n->id+1 << ")";
			//}
			
			if(fabs(Sgrains[i].x-MIN)<TOL2) {F[0].push_back(i+1); rigidF[0].insert(compid+1);}
			if(fabs(Sgrains[i].y-MIN)<TOL2) {F[1].push_back(i+1); rigidF[1].insert(compid+1);}
			if(fabs(Sgrains[i].z-MIN)<TOL2) {F[2].push_back(i+1); rigidF[2].insert(compid+1);}
			
			if(fabs(Sgrains[i].x-MAX)<TOL2) {L[0].push_back(i+1); rigidL[0].insert(compid+1);}
			if(fabs(Sgrains[i].y-MAX)<TOL2) {L[1].push_back(i+1); rigidL[1].insert(compid+1);}
			if(fabs(Sgrains[i].z-MAX)<TOL2) {L[2].push_back(i+1); rigidL[2].insert(compid+1);}
		}
		
		outLfile << "\nBonds\n\n";
		count_bonds = 0;
		for(int i=0;i<Sedges.size();i++)
	    {
			//if(component[Sedges[i].m->id]==component[Sedges[i].n->id])
			if(Sedges[i].active)
			{
				outLfile << count_bonds+1 << " 1 " << Sedges[i].m->id+1 << " " << Sedges[i].n->id+1 << std::endl;
				count_bonds++;
			}
			
			//outLfile << i+1 << " 1 " << Sedges[i].m->id+1 << " " << Sedges[i].n->id+1 << std::endl;
		}
		
		std::vector<std::pair<int, int>> S1 = { { 0, 1 }, { 1, 3 }, { 3, 2 }, { 2, 0 }, { 4, 5 }, { 5, 7 }, { 7, 6 }, { 6, 4 }, { 8, 9 }, { 9, 11 }, { 11, 10 }, { 10, 8 }};
		std::vector<std::pair<int, int>> S2 = { { 6, 10 }, { 10, 7 }, { 7, 11 }, { 11, 6 }, { 4, 8 }, { 8, 5 }, { 5, 9 }, { 9, 4 },
												{ 2, 9 }, { 9, 3 }, { 3, 11 }, { 11, 2 }, { 0, 8 }, { 8, 1 }, { 1, 10 }, { 10, 0 },
												{ 1, 5 }, { 5, 3 }, { 3, 7 }, { 7, 1 }, { 0, 4 }, { 4, 2 }, { 2, 6 }, { 6, 0 } };
		
		outLfile << "\nAngles\n\n";
		
		//std::cout << std::endl << "Map";
		int count_angles = 0;
		for(int i=0;i<Sgrains.size();i++)
	    {
			std::vector<int> array0(12,-1);
			for(int j=0;j<Sgrains[i].sfaces.size();j++)
			{
				int iter = Sgrains[i].sfaces[j];
				int other_id = (Sedges[iter].m->id == i) ? Sedges[iter].n->id : Sedges[iter].m->id;
				
				//if(component[other_id]==component[i])
				if(Sedges[iter].active)
				{
					double d[3] = {Sgrains[other_id].x - Sgrains[i].x, Sgrains[other_id].y - Sgrains[i].y, Sgrains[other_id].z - Sgrains[i].z};
					double dd = sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
					for(int i2=0;i2<3;i2++) d[i2] = std::round(d[i2]*sqrt(2)/dd);
					
					double face_index;
					
					if(fabs(d[0])<TOL2) face_index = (d[1]+1) + 0.5*(d[2]+1);
					if(fabs(d[1])<TOL2) face_index = 4 + (d[0]+1) + 0.5*(d[2]+1);
					if(fabs(d[2])<TOL2) face_index = 8 + (d[0]+1) + 0.5*(d[1]+1);
					
					std::cout << std::endl << d[0] << " " << d[1] << " " << d[2] << " " << face_index;
					int findex = (int) face_index;
					
					array0[findex] = other_id;
				}
			}
			
			for( std::vector<std::pair<int, int>>::iterator iter0 = S1.begin(); iter0 < S1.end(); iter0++ )
			{
				int id1 = iter0->first;
				int id2 = iter0->second;
				
				if( (array0[id1]!=-1) && (array0[id2]!=-1) )
				{
					outLfile << count_angles+1 << " 1 " << array0[id1]+1 << " " << i+1 << " " << array0[id2]+1 << std::endl;
					count_angles++;
				}
			}
			
			for( std::vector<std::pair<int, int>>::iterator iter0 = S2.begin(); iter0 < S2.end(); iter0++ )
			{
				int id1 = iter0->first;
				int id2 = iter0->second;
				
				if( (array0[id1]!=-1) && (array0[id2]!=-1) )
				{
					outLfile << count_angles+1 << " 2 " << array0[id1]+1 << " " << i+1 << " " << array0[id2]+1 << std::endl;
					count_angles++;
				}
			}
			
		}
		
		outLfile.close();
		std::ifstream infile1("ABA.temp");
		
		outLfile2 << "LAMMPS data file\n\n\t" << Sgrains.size() << " atoms\n\t"<< count_bonds <<" bonds\n\t" << count_angles << " angles\n\t0 dihedrals\n\t0 impropers\n\n";
	  outLfile2 << "\t1 atom types\n\t1 bond types\n\t2 angle types\n\t0 dihedral types\n\t0 improper types\n\n";
	  outLfile2 << "# Cell: orthorhombic\n\t"<< -sizeN*0.5000 << "\t" << sizeN*1.5000 << " xlo xhi\n\t";
	  outLfile2 << -sizeN*0.5000 << "\t" << sizeN*1.5000 << " ylo yhi\n\t";
	  outLfile2 << -sizeN*0.5000 << "\t" << sizeN*1.5000 << " zlo zhi\n\n";
	  //outLfile << "# Cell: orthorhombic\n\t-5.00000000\t10.00000000 xlo xhi\n\t"<< MIN-(0.25/sizeN) << "\t" << MAX+(0.25/sizeN) << " ylo yhi\n\t" << MIN-(0.25/sizeN) << "\t" << MAX+(0.25/sizeN) << " zlo zhi\n\n";
	  outLfile2 << "Masses\n\n\t1 27.00 \n\n";
		
		outLfile2 << "Bond Coeffs\n\n 1 18.71123 2.82842712 \n\n";
		outLfile2 << "Angle Coeffs\n\n 1 8.63595 90.0 # 90 deg.\n 2 2.87865 60.0 # 60 deg.\n";
				
		while(std::getline(infile1, line))
		{
			outLfile2 << line << std::endl;;
		}
		
		infile1.close();
		
		//std::cout << std::endl << "End Map";
		
		for(int i=0;i<3;i++)
		{
			if(F[i].size()>0)
			{
				int countF=0;
				//outLfileG << "group F" << i << " id"; 
				for(std::vector<int>::iterator it = F[i].begin(); it < F[i].end(); ++it,++countF)
				{
					if(countF%32==0) outLfileG << "\ngroup F" << i << " id"; 
					outLfileG << " " << *it;
					
				}
			}
			
			if(L[i].size()>0)
			{
				int countL=0;
				//outLfileG << "group L" << i << " id"; 
				for(std::vector<int>::iterator it = L[i].begin(); it < L[i].end(); ++it,++countL)
				{
					if(countL%32==0) outLfileG << "\ngroup L" << i << " id"; 
					outLfileG << " " << *it;
					
				}
			}
			
			if(rigidF[i].size()>0)
			{
				int countF=0;
				//outLfileG << "group F" << i << " id"; 
				for(std::set<int>::iterator it = rigidF[i].begin(); it != rigidF[i].end(); ++it,++countF)
				{
					if(countF%32==0) outLfileG << "\ngroup rigidF" << i << " molecule"; 
					outLfileG << " " << *it;
					
				}
			}
			
			if(rigidL[i].size()>0)
			{
				int countL=0;
				//outLfileG << "group F" << i << " id"; 
				for(std::set<int>::iterator it = rigidL[i].begin(); it != rigidL[i].end(); ++it,++countL)
				{
					if(countL%32==0) outLfileG << "\ngroup rigidL" << i << " molecule"; 
					outLfileG << " " << *it;
					
				}
			}
		}
		
		std::vector<int> rigid_objects;
		int max=0,comp_max=0;
		for(int i=0;i<comp.size();i++)
		{
			//std::cout << std::endl << " Component " << i << " has " << comp[i].contents.size() << " grains :";
			
			if( comp[i].contents.size() >= 2)
			{
				rigid_objects.push_back(i+1);
				
				//max = comp[i].contents.size();
				//comp_max = i;
				//std::cout << "\n Hello\n";
			}
			/*
			for(int j=0;j<comp[i].contents.size();j++)
			{
				std::cout << " " << comp[i].contents[j];
			}*/
		}
		
		if(rigid_objects.size()>0)
		{
			int countL=0;
			//outLfileG << "group L" << i << " id"; 
			for(std::vector<int>::iterator it = rigid_objects.begin(); it != rigid_objects.end(); ++it,++countL)
			{
				if(countL%32==0) outLfileG << "\ngroup rigid_objects molecule"; 
				outLfileG << " " << *it;
				
			}
		}
		
		//std::cout << std::endl << "Largest Component " << comp_max << " which has " << max << " grains\n\n";
		
		///////////////////////// END WRITING LAMMPS FILES (data.lmp & groups.dat) ////////////////////////////
		
		/*
		int max=0,comp_max=0;
		for(int i=0;i<comp.size();i++)
		{
			std::cout << std::endl << " Component " << i << " has " << comp[i].contents.size() << " grains :";
			if( comp[i].contents.size() >= max)
			{
				max = comp[i].contents.size();
				comp_max = i;
				//std::cout << "\n Hello\n";
			}
			for(int j=0;j<comp[i].contents.size();j++)
			{
				std::cout << " " << comp[i].contents[j];
			}
		}
		
		std::cout << std::endl << "Largest Component " << comp_max << " which has " << max << " grains\n\n";
		*/
		
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
		
		std::stringstream bufferc;
		bufferc << "solids/Assembly_all.stl.cracks";
		std::ofstream outfile2(bufferc.str().c_str());
		
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
											
						if(  (cid2!=i) || ( (cid2==i)&&(curr_grain<grain2) )  ) // this ensures only faces across different components are treated AND internal faces are treated exactly only once for the smaller internal component
						{
							double dir[3] = {Sgrains[grain2].x-Sgrains[curr_grain].x,Sgrains[grain2].y-Sgrains[curr_grain].y,Sgrains[grain2].z-Sgrains[curr_grain].z};
							double dirmag = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
							double dirn[3] = {dir[0]/dirmag,dir[1]/dirmag,dir[2]/dirmag};
							
							double u[3] = {Sedges[curr_edge].vertices[3]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[4]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[5]-Sedges[curr_edge].vertices[2]}; 
							double v[3] = {Sedges[curr_edge].vertices[6]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[7]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[8]-Sedges[curr_edge].vertices[2]}; 
							double uxv[3] = {u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
							double nmag = sqrt(uxv[0]*uxv[0] + uxv[1]*uxv[1] + uxv[2]*uxv[2]);
							double n[3] = {uxv[0]/nmag,uxv[1]/nmag,uxv[2]/nmag};
							
							double dot = dirn[0]*n[0]+dirn[1]*n[1]+dirn[2]*n[2]; // dot=+1 or -1
							//outfile << std::endl << "Dot " << dot;
							
							int order[3] = {0,1,2};
							if(dot<0) {order[1]=2; order[2]=1;}
							
							if (Sedges[curr_edge].vertices.size()==9)
							{
								
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
							else
							{
								
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
								
								order[1] += 1; order[2] += 1;
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
						
						int order[3] = {0,1,2};
						if(dot<0) {order[1]=2; order[2]=1;}
						
						if (Bedges[curr_edge].vertices.size()==9)
						{
							
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
						else
						{
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
							
							order[1] += 1; order[2] += 1;
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
						////
					
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
								
					//looping over internal faces
					for(int k=0;k<Sgrains[curr_grain].sfaces.size();k++)
					{
						int curr_edge = Sgrains[curr_grain].sfaces[k];
						
						if(Sedges[curr_edge].m->id == curr_grain)
							grain2 = Sedges[curr_edge].n->id;
						else if(Sedges[curr_edge].n->id == curr_grain)
							grain2 = Sedges[curr_edge].m->id;
						
						int cid2 = component[grain2];
						
						//if(cid2 != i) * removed because in the inclusion all internal faces must be present so that crack faces in the .cracks can be matched with mesh faces that have been previously defined
						if(  (cid2!=i) || ( (cid2==i)&&(curr_grain<grain2) )  ) // this ensures only faces across different components are treated AND internal faces are treated exactly only once
						{
							double dir[3] = {Sgrains[grain2].x-Sgrains[curr_grain].x,Sgrains[grain2].y-Sgrains[curr_grain].y,Sgrains[grain2].z-Sgrains[curr_grain].z};
							double dirmag = sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
							double dirn[3] = {dir[0]/dirmag,dir[1]/dirmag,dir[2]/dirmag};
							
							double u[3] = {Sedges[curr_edge].vertices[3]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[4]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[5]-Sedges[curr_edge].vertices[2]}; 
							double v[3] = {Sedges[curr_edge].vertices[6]-Sedges[curr_edge].vertices[0],Sedges[curr_edge].vertices[7]-Sedges[curr_edge].vertices[1],Sedges[curr_edge].vertices[8]-Sedges[curr_edge].vertices[2]}; 
							double uxv[3] = {u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]};
							double nmag = sqrt(uxv[0]*uxv[0] + uxv[1]*uxv[1] + uxv[2]*uxv[2]);
							double n[3] = {uxv[0]/nmag,uxv[1]/nmag,uxv[2]/nmag};
							
							double dot = dirn[0]*n[0]+dirn[1]*n[1]+dirn[2]*n[2]; // dot=+1 or -1
							//outfile << std::endl << "Dot " << dot;
							
							int order[3] = {0,1,2};
							if(dot<0) {order[1]=2; order[2]=1;}
							
							//if(cid2!=i) * removed because in the inclusion all internal faces must be present so that crack faces can match with mesh faces that have been previously defined
							outfile << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
							outfile_solid << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
							
							//outfile1 << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
																				
							for(int i1=0; i1<3; i1++)
							{
								int vloc = vertex_names.at(   Sedges[curr_edge].vertex_strs[(order[i1])]   );
								
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
								//outfile << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
								outfile << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
								outfile_solid << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
																
								if((Sedges[curr_edge].active==0)&&(cid2==i)) outfile2 << "(" << Svertices[vloc].x-TOL*toln[0] << "," << Svertices[vloc].y-TOL*toln[1] << "," << Svertices[vloc].z-TOL*toln[2] << ") ";
								//outfile1 << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2] << " " << Svertices[vloc].id << "Vertex " << Sedges[curr_edge].vertex_strs[(order[i1])] << " " << Svertices[vloc].id << " " << Svertices[vloc].x << " " << Svertices[vloc].y << " " << Svertices[vloc].z << " (" << Svertices[vloc].grains.size() << ") ";
								//std::cout << std::endl << Svertices[vloc].id << "      Fvertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
							}
							if((Sedges[curr_edge].active==0)&&(cid2==i)) outfile2 << std::endl;
							
							//if(cid2!=i) 
							outfile << std::endl << "    endloop" << std::endl << "  endfacet";
							outfile_solid << std::endl << "    endloop" << std::endl << "  endfacet";
							//outfile1 << std::endl << "    endloop" << std::endl << "  endfacet";
							
							if (Sedges[curr_edge].vertices.size()==12) // EXTRA triangular FACE if rhombus
							{
								order[1] += 1; order[2] += 1;
								
								//if(cid2!=i) 
								outfile << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
								outfile_solid << std::endl << "  facet normal  " << dirn[0] << "  " << dirn[1] << "  " << dirn[2] << std::endl << "    outer loop";
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
									//outfile << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
									outfile << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
									outfile_solid << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
									if((Sedges[curr_edge].active==0)&&(cid2==i)) outfile2 << "(" << Svertices[vloc].x-TOL*toln[0] << "," << Svertices[vloc].y-TOL*toln[1] << "," << Svertices[vloc].z-TOL*toln[2] << ") ";
									//outfile1 << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2] << " " << Svertices[vloc].id << "Vertex " << Sedges[curr_edge].vertex_strs[(order[i1])] << " " << Svertices[vloc].id << " " << Svertices[vloc].x << " " << Svertices[vloc].y << " " << Svertices[vloc].z << " (" << Svertices[vloc].grains.size() << ") ";
									//std::cout << std::endl << Svertices[vloc].id << "      Fvertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
								}
								
								//if(cid2!=i) 
								outfile << std::endl << "    endloop" << std::endl << "  endfacet";
								outfile_solid << std::endl << "    endloop" << std::endl << "  endfacet";
								
								if((Sedges[curr_edge].active==0)&&(cid2==i)) outfile2 << std::endl;
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
						
						int order[3] = {0,1,2};
						if(dot<0) {order[1]=2; order[2]=1;}
						
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
							
							//outfile << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
							outfile << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
							outfile_solid << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
							//outfile1 << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2] << " " << Svertices[vloc].id << "Vertex " << Bedges[curr_edge].vertex_strs[(order[i1])] << " " << Svertices[vloc].id << " " << Svertices[vloc].x << " " << Svertices[vloc].y << " " << Svertices[vloc].z << " (" << Svertices[vloc].grains.size() << ") ";
						}
						
						
						outfile << std::endl << "    endloop" << std::endl << "  endfacet";
						outfile_solid << std::endl << "    endloop" << std::endl << "  endfacet";
						//outfile1 << std::endl << "    endloop" << std::endl << "  endfacet";
						
						if (Bedges[curr_edge].vertices.size()==12) // EXTRA triangular FACE if rhombus
						{
							order[1] += 1; order[2] += 1;
							
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
								
								//outfile << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2];
								outfile << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
								outfile_solid << std::endl << "      vertex  " << std::round((Svertices[vloc].x-TOL*toln[0]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].y-TOL*toln[1]-MIN)*sizeN-2) << " " << std::round((Svertices[vloc].z-TOL*toln[2]-MIN)*sizeN-2);
								//outfile1 << std::endl << "      vertex  " << Svertices[vloc].x-TOL*toln[0] << " " << Svertices[vloc].y-TOL*toln[1] << " " << Svertices[vloc].z-TOL*toln[2] << " " << Svertices[vloc].id << "Vertex " << Bedges[curr_edge].vertex_strs[(order[i1])] << " " << Svertices[vloc].id << " " << Svertices[vloc].x << " " << Svertices[vloc].y << " " << Svertices[vloc].z << " (" << Svertices[vloc].grains.size() << ") ";
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
				/*
				for(std::map<std::string, int>::const_iterator it = vertex_names.begin(); it != vertex_names.end(); ++it)
				{
					std::cout << it->first << " " << it->second << " " << Svertices[it->second].x << " " << Svertices[it->second].y << " " << Svertices[it->second].z << "\n";
				}*/
				
				//outfile << std::endl << "endsolid" << std::endl;
				//outfile1 << std::endl << "endsolid" << std::endl;
			
		}
		
		outfile.close();
		outfile2.close();
			
		/*
		for(std::map<std::string, int>::const_iterator it = vertex_names.begin(); it != vertex_names.end(); ++it)
		{
			std::cout << it->first << " " << it->second << " " << Svertices[it->second].x << " " << Svertices[it->second].y << " " << Svertices[it->second].z << "\n";
		}
		
		for(int i=0;i<Svertices.size();i++)
		{
				outfile << Svertices[i].x << " " << Svertices[i].y << " " << Svertices[i].z << "\n";
		}
		
		for(int i=0;i<faces.size();i++)
		{
				outfile1 << faces[i];
				if(i%3==2) outfile1 << std::endl;
				else outfile1 << " ";
		}
		*/
		//outfile1.close();
		//outfile.close();
	    
	}
	
	std::cout << std::endl;
	
	
	
	infile.close();
	//outfile.close();
	
	int stop_s=clock();
	std::cout << "Total time: " << (stop_s-start_s)/double(CLOCKS_PER_SEC) << std::endl;
	
	return 0;
	
}

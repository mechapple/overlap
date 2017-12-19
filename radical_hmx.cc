// Radical Voronoi tessellation example code
//
// Author   : Chris H. Rycroft (LBL / UC Berkeley)
// Email    : chr@alum.mit.edu
// Date     : August 30th 2011

#include "voro++.hh"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace voro;

// Set up constants for the container geometry
//const double xlo=0,xhi=128.02344000;
//const double ylo=0,yhi=123.05028000;
//const double zlo=0,zhi=7.5143014286;
//const double xy=0.0,xz=-1.1242628571,yz=0.0;

// Set up the number of blocks that the container is divided
// into.
const int n_x=6,n_y=6,n_z=5;

int main() 
{

	std::ifstream infileb("bounds.dat");
	
	std::string line,input;
	std::stringstream sstream;
	
	double xlo=0.0,xhi=1.0;
	double ylo=0.0,yhi=1.0;
	double zlo=0.0,zhi=1.0;
	double xy=0.0,xz=0.0,yz=0.0;
	int countb=0;
	
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
		if(countb==1) { xlo = ::atof(entries[0].c_str()); xhi = ::atof(entries[1].c_str());}
		if(countb==2) { ylo = ::atof(entries[0].c_str()); yhi = ::atof(entries[1].c_str());}
		if(countb==3) { zlo = ::atof(entries[0].c_str()); zhi = ::atof(entries[1].c_str());}
		if(countb==4) { xy = ::atof(entries[0].c_str()); xz = ::atof(entries[1].c_str()); yz = ::atof(entries[2].c_str());}
	}
	
	printf("\n %g %g %g %g %g %g %g %g %g\n",xlo,xhi,ylo,yhi,zlo,zhi,xy,xz,yz);
	
	// Create a container with the geometry given above, and make it
	// non-periodic in each of the three coordinates. Allocate space for
	// eight particles within each computational block. Import
	// the monodisperse test packing and output the Voronoi
	// tessellation in gnuplot and POV-Ray formats.
	//container con(xlo,xhi,ylo,yhi,zlo,zhi,n_x,n_y,n_z,
			//false,false,false,8);
	//con.import("pack_six_cube");
	//con.draw_cells_gnuplot("pack_six_cube.gnu");
	//con.draw_cells_pov("pack_six_cube_v.pov");
	//con.draw_particles_pov("pack_six_cube_p.pov");

	// Create a container for polydisperse particles using the same
	// geometry as above. Import the polydisperse test packing and
	// output the Voronoi radical tessellation in gnuplot and POV-Ray
	// formats.
	
	//container_poly conp(xlo,xhi,ylo,yhi,zlo,zhi,n_x,n_y,n_z,false,false,false,8);
	
	container_periodic_poly conp(xhi-xlo,xy,yhi-ylo,xz,yz,zhi-zlo,n_x,n_y,n_z,8);
	
	conp.import("hmx.dat");
	conp.draw_cells_gnuplot("hmx.gnu");
	conp.draw_cells_pov("hmx_v.pov");
	conp.draw_particles_pov("hmx_p.pov");
	double vvol=conp.sum_cell_volumes();
	double cvol=(xhi-xlo)*(yhi-ylo)*(zhi-zlo);
	
	std::cout << "\n Volumes " << vvol << " " << cvol << " packing fraction " << vvol/cvol << "\n";
	//%i The particle ID number.
	//%q The position vector of the particle, short for “%x %y %z”.
	//%s The number of faces of the Voronoi cell.
	//%n A list of the neighboring particle or wall IDs corresponding to each face.
	//%f A list of areas of each face.
	//%t A list of bracketed sequences of vertices that make up each face.
	//%w The number of vertices in the Voronoi cell.
	//%P A list of the vertices of the Voronoi cell in the format (x,y,z), relative to the global coordinate system
	conp.print_custom("%i %q %s %n %f %t %w %P","voro_cell.dat");
	
	infileb.close();
}

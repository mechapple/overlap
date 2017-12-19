//#include "mex.h"
//g++ -frounding-math -std=c++11 get_intersection_mesh.cpp -I/opt/libigl/include 
//			-I/opt/eigen -lboost_thread -lCGAL -lmpfr -lgmp -lboost_system
// ./a.out <ngridx> <ngridy>

//#include <igl/viewer/Viewer.h>
#include <Eigen/Core>
//#include <igl/matlab/prepare_lhs.h>
//#include <igl/matlab/parse_rhs.h>
#include <igl/copyleft/cgal/mesh_boolean.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#define MAXBUFSIZE  ((int) 1e6)

int sgn(double x)  //1 if x>=0 and -1 if x<0
{
	int val = (x >= 0) - (x < 0);
	return val;
}

struct point{
	double x,y,z;
	point(double a, double b, double c) : x(a), y(b), z(c) {}
};

struct boundary {
	double xlo,xhi;
	double ylo,yhi;
	double zlo,zhi;
	double xy,xz,yz;
};

double_t SignedVolumeOfTriangle(point p1, point p2, point p3) {
    double v321 = p3.x*p2.y*p1.z;
    double v231 = p2.x*p3.y*p1.z;
    double v312 = p3.x*p1.y*p2.z;
    double v132 = p1.x*p3.y*p2.z;
    double v213 = p2.x*p1.y*p3.z;
    double v123 = p1.x*p2.y*p3.z;
    return (1.0/6.0)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

Eigen::MatrixXd readMatrixd(const char *filename)
{
	int cols = 0, rows = 0;
	double buff[MAXBUFSIZE];
	
	// Read numbers from file into buffer.
	std::ifstream infile;
	infile.open(filename);
	while (! infile.eof())
	{
		std::string line;
		getline(infile, line);
	
		int temp_cols = 0;
		std::stringstream stream(line);
		while(! stream.eof())
			stream >> buff[cols*rows+temp_cols++];
	
		if (temp_cols == 0)
			continue;
	
		if (cols == 0)
			cols = temp_cols;
	
		rows++;
	}
	
	infile.close();
	
	rows--;
	
	// Populate matrix with numbers.
	Eigen::MatrixXd result(rows,cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result(i,j) = buff[ cols*i+j ];
	
	return result;
};

Eigen::MatrixXi readMatrixi(const char *filename)
{
	int cols = 0, rows = 0;
	double buff[MAXBUFSIZE];
	
	// Read numbers from file into buffer.
	std::ifstream infile;
	infile.open(filename);
	while (! infile.eof())
	{
		std::string line;
		getline(infile, line);
	
		int temp_cols = 0;
		std::stringstream stream(line);
		while(! stream.eof())
			stream >> buff[cols*rows+temp_cols++];
	
		if (temp_cols == 0)
			continue;
	
		if (cols == 0)
			cols = temp_cols;
	
		rows++;
	}
	
	infile.close();
	
	rows--;
	
	// Populate matrix with numbers.
	Eigen::MatrixXi result(rows,cols);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result(i,j) = buff[ cols*i+j ];
	
	return result;
};

int main(int argc, char *argv[])
{
	Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
	std::string sep = "\n----------------------------------------\n";
	
	//******** BOUNDARY INFO ***********//
	std::string line,input;
	std::stringstream sstream;
	std::ifstream infileb("bounds.cell");
	boundary bounds;
	
	double dx = 1.0/atof(argv[1]), dy = 1.0/atof(argv[2]);
	double TOL = 1e-6;
	
	//std::cout << "\nIncrements " << dx << " " << dy << "\n";
	//double xlo=0.0,xhi=1.0;
	//double ylo=0.0,yhi=1.0;
	//double zlo=0.0,zhi=1.0;
	//double xy=0.0,xz=0.0,yz=0.0;
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
		if(countb==1) { bounds.xlo = ::atof(entries[0].c_str()); bounds.xhi = ::atof(entries[1].c_str());}
		if(countb==2) { bounds.ylo = ::atof(entries[0].c_str()); bounds.yhi = ::atof(entries[1].c_str());}
		if(countb==3) { bounds.zlo = ::atof(entries[0].c_str()); bounds.zhi = ::atof(entries[1].c_str());}
		if(countb==4) { bounds.xy = ::atof(entries[0].c_str()); bounds.xz = ::atof(entries[1].c_str()); bounds.yz = ::atof(entries[2].c_str());}
	}
	
	//printf("\n %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g %4.16g\n",bounds.xlo,bounds.xhi,bounds.ylo,bounds.yhi,bounds.zlo,bounds.zhi,bounds.xy,bounds.xz,bounds.yz);
	
	//double xprd=(xhi-xlo), yprd=(yhi-ylo), zprd=(zhi-zlo);
	infileb.close();
	//**********************************//
	
	
    Eigen::MatrixXd V1 = readMatrixd("V1.dat");
    Eigen::MatrixXi F1 = readMatrixi("F1.dat");
    
    //std::cout << V1.format(HeavyFmt) << sep;
    //std::cout << F1 << sep;
    
    Eigen::MatrixXd V2o = readMatrixd("V2.dat");
    Eigen::MatrixXi F2 = readMatrixi("F2.dat");
    
    //std::cout << V2.format(HeavyFmt) << sep;
    //std::cout << F2 << sep;
    
    for(double ix=-0.5;ix<=(0.5+TOL);ix+=dx)
		for(double iy=-0.5;iy<=(0.5+TOL);iy+=dy)
		{
    
			Eigen::MatrixXd disp = Eigen::MatrixXd::Zero(V2o.rows(),V2o.cols());
		    for(size_t i=0; i<disp.rows(); i++)
		    {
				disp(i,0) = (ix*(bounds.xhi-bounds.xlo)+iy*(bounds.xy)); //*sgn(-bounds.xz) ;
				disp(i,1) = iy*(bounds.yhi-bounds.ylo); //*sgn(-bounds.yz);
			}
		    
		    Eigen::MatrixXd V2=V2o+disp;
		    		    
		    Eigen::MatrixXd V_i;
			Eigen::MatrixXi F_i;
		  
			igl::copyleft::cgal::mesh_boolean(V1,F1,V2,F2,igl::MESH_BOOLEAN_TYPE_INTERSECT,V_i,F_i);
		  
			double sum_vol2=0,sum_voli=0;
		    for(size_t i=0; i<F2.rows(); i++)
		    {
				point p1( V2(F2(i,0),0), V2(F2(i,0),1), V2(F2(i,0),2) );
				point p2( V2(F2(i,1),0), V2(F2(i,1),1), V2(F2(i,1),2) );
				point p3( V2(F2(i,2),0), V2(F2(i,2),1), V2(F2(i,2),2) );
				
				sum_vol2 += SignedVolumeOfTriangle(p1,p2,p3);
			}
		    
			for(size_t i=0; i<F_i.rows(); i++)
		    {
				point p1( V_i(F_i(i,0),0), V_i(F_i(i,0),1), V_i(F_i(i,0),2) );
				point p2( V_i(F_i(i,1),0), V_i(F_i(i,1),1), V_i(F_i(i,1),2) );
				point p3( V_i(F_i(i,2),0), V_i(F_i(i,2),1), V_i(F_i(i,2),2) );
				
				sum_voli += SignedVolumeOfTriangle(p1,p2,p3);
			}
		  
			//std::cout << "\nVol2 = " << sum_vol2;
			//std::cout << "\nVoli = " << sum_voli << "\n";
			
			printf("%g %g %g\n",ix,iy,sum_voli);
		}
	
	//printf("\n");
	
	//std::cout << V_i.format(HeavyFmt) << sep;
	//std::cout << F_i.format(HeavyFmt) << sep;
  
	//igl::viewer::Viewer viewer;
	//viewer.data.set_mesh(V, F);
	//viewer.data.set_face_based(true);
	//viewer.launch();
	
	return 0;
}

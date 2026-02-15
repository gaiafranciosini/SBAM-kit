 /*
 *  
 *  Fred project
 *
 *  
 *  A. Schiavi (2020)
 */
#include <string>
#include <iostream>
using namespace std;

class Voxels { // wrapper class
public:	
	string dtype;
	char *data; // opaque pointer to voxel values 
	int nn[3] ; // dims
	float x0[3]; // "lowest corner" in global coords
	float offset[3]; // "centroid of first voxel" in global coords 
	float L[3],hs[3]; // extent and spacing
	float left[3],up[3],front[3]; // orthonormal basis wrt offset point, i.e. the origin
	bool good;
	const size_t &N; // read-only number of voxels
private:
	size_t _N; // internal number of voxels
	void hs2L(){for(int i=0;i<3;i++) L[i]=nn[i]*hs[i];}
	void L2hs(){for(int i=0;i<3;i++) hs[i]=L[i]/nn[i];}
public:
	Voxels();
	Voxels(string dtype,int nn[3],float hs[3], float x0[3]);

	~Voxels();

	void reset(string dtype,int nn[3],float hs[3], float x0[3]);


	int read(string path);
	int write(string path,bool localData=true);

	void info(ostream &os = cout);

	size_t dataSize();

	void clone(Voxels &B);

	void resize(int nnnew[3]);
	void resize(int nx,int ny,int nz);

	void resize_with_spacing(float hs[3]);
	void resize_with_spacing(float hx,float hy,float hz);

	void eye();

	void realloc(string dtype,int nn[3]);

	void setSpacing(float hs[3]);
	void setLowestCorner(float x0[3]);
	void updateOffset();
	void setBasis(float left[3],float up[3],float front[3]);

};


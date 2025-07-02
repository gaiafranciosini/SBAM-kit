#include <iostream>
#include <cmath>
#include "types.h"


using namespace std;

#pragma once

template <typename T>
class idx3d {
public:
    union {
		T v[3];
		struct { T x,y,z;}; // CART
		struct { T Alpha,R,Z;}; // CYL
		struct { T phi,theta,rho;}; // SPHER
	};

	idx3d(){};

	idx3d(T a,T b, T c){
 		v[0]=a;v[1]=b;v[2]=c;
 	};
	
	idx3d(const idx3d &a){
 		v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];
 	};
	
	idx3d(T n){
 		v[0]=v[1]=v[2]=n;
 	}
	
 	idx3d(T *a){
 		v[0]=a[0];v[1]=a[1];v[2]=a[2];
 	}
 	
	inline T & operator[](int j){return v[j];};
	inline const T & operator[](int j) const {return v[j];};
	
	inline void set(T x,T y,T z){v[0]=x;v[1]=y;v[2]=z;}
	
	
	inline idx3d& operator=(const idx3d &a){v[0]=a.v[0];v[1]=a.v[1];v[2]=a.v[2];return *this;}
    
    inline idx3d& operator+=(T a){v[0]+=a;v[1]+=a;v[2]+=a;return *this;};
    inline idx3d& operator-=(T a){v[0]-=a;v[1]-=a;v[2]-=a;return *this;};
    inline idx3d& operator*=(T a){v[0]*=a;v[1]*=a;v[2]*=a;return *this;};
    inline idx3d& operator/=(T a){v[0]/=a;v[1]/=a;v[2]/=a;return *this;};
    inline idx3d&  operator=(T a){v[0] =a;v[1] =a;v[2] =a;return *this;};
    
    inline idx3d& operator++(){v[0]++;v[1]++;v[2]++;return *this;};
    inline idx3d& operator--(){v[0]--;v[1]--;v[2]--;return *this;};
    inline idx3d& operator++(int ){v[0]++;v[1]++;v[2]++;return *this;};
    inline idx3d& operator--(int ){v[0]--;v[1]--;v[2]--;return *this;};
    
    inline idx3d& operator+=(const idx3d &a){v[0]+=a.v[0];v[1]+=a.v[1];v[2]+=a.v[2];return *this;}
    inline idx3d& operator-=(const idx3d &a){v[0]-=a.v[0];v[1]-=a.v[1];v[2]-=a.v[2];return *this;}
    inline idx3d& operator*=(const idx3d &a){v[0]*=a.v[0];v[1]*=a.v[1];v[2]*=a.v[2];return *this;}
    inline idx3d& operator/=(const idx3d &a){v[0]/=a.v[0];v[1]/=a.v[1];v[2]/=a.v[2];return *this;}
    
    inline idx3d operator-() const {return idx3d(-v[0],-v[1],-v[2]);}
    
    friend idx3d operator+(idx3d lhs,       // passing first arg by value helps optimize chained a+b+c
                    const idx3d& rhs) // alternatively, both parameters may be const references.
  	{
    	 return lhs += rhs; // reuse compound assignment and return the result by value
  	}
    friend idx3d operator-(idx3d lhs,const idx3d& rhs){return lhs -= rhs;}
    friend idx3d operator*(idx3d lhs,const idx3d& rhs){return lhs *= rhs;}
    friend idx3d operator/(idx3d lhs,const idx3d& rhs){return lhs /= rhs;}
	
};

template <typename T>
ostream& operator<<(ostream &os,const idx3d<T> &v){
	os<<'['<<v[0]<<','<<v[1]<<','<<v[2]<<']';
	return os;
}

template <typename T>
idx3d<T> min(const idx3d<T> &a, const idx3d<T> &b){
  idx3d<T> x;
  for (int j=0;j<3;j++) x[j]=a[j]<b[j]? a[j]: b[j];
  return x;
}

template <typename T>
idx3d<T> max(const idx3d<T> &a, const idx3d<T> &b){
  idx3d<T> x;
  for (int j=0;j<3;j++) x[j]=a[j]>b[j]? a[j]: b[j];
  return x;
}

template <typename T>
bool operator==(const idx3d<T> &a,const idx3d<T> &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return false;
	return true;
}

template <typename T>
bool operator!=(const idx3d<T> &a,const idx3d<T> &b){
    if(a.x != b.x || a.y != b.y || a.z != b.z ) return true;
	return false;
}


typedef idx3d<int32> i3d;
typedef idx3d<int16> i3ds;
typedef idx3d<uint16> ui3ds;

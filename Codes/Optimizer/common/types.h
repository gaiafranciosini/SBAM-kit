/*
 *  types.h
 *  Fred
 *
 *  Created by Angelo Schiavi on 18/01/08.
 *  Copyright 2008 A. Schiavi. All rights reserved.
 *
 */

#pragma once

typedef signed char           int8;       // 8 bit signed
typedef unsigned char         uint8;      // 8 bit unsigned
typedef short                 int16;      // 16 bit signed
typedef unsigned short        uint16;     // 16 bit unsigned
typedef int                   int32;      // 32 bit signed
typedef unsigned int          uint32;     // 32 bit unsigned
typedef long long             int64;      // 64 bit signed
typedef unsigned long long    uint64;     // 64 bit unsigned

typedef uint32 uint;

typedef float                 float32;     // 32 bit floating point 
typedef double                float64;     // 64 bit floating point

// new particle types must be added here and in Particle.cpp
const int32 rayTypeSize = 4;  // NB: must be equal to rayType size!!!
enum rayType {
	GHOST=0,GEORAY,H1,C12
};

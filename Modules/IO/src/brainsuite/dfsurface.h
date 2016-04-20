// Copyright (c) 2010 David W. Shattuck, PhD
// $Id: dfsurface.h 154 2010-08-24 22:15:01Z shattuck $
#ifndef DFSurface_H
#define DFSurface_H

/**
 \brief  BrainSuite Surface Format file reader classes
 \author David Shattuck 
 \date   $Date: 2010-08-24 15:15:01 -0700 (Tue, 24 Aug 2010) $
 This file contains C++ code to read and write the BrainSuite surface file format (dfs).
 
 A DFS file contains the following layout:

 [ DFS Header           ] <br>
 [ metadata             ] embedded XML (not used) <br>
 [ subject data         ] embedded XML (not used) <br>
 [ triangle data        ] a list of triangles, specified by an (a,b,c) triple of 32-bit integer numbers<br>
 [ vertices             ] a list of vertices, specified by an (x,y,z) triple of 32-bit floating point numbers<br>
 [ vertex normals       ]	vertex normal data (0 if not in file) (nx,ny,nz) in 32-bit floating point<br>
 [ vertex UV coordinates] surface parameterization data (u,v) in 32-bit floating point<br>
 [ vertex colors        ]	per vertex color data in (r,g,b) format in 32-bit floating point ([0-1])<br>
 [ vertex labels        ] per vertex label index specified in 16-bit integer format<br>
 [ vertex attributes    ]	vertex attributes (float32 array of length NV)<br>

 The header portion contains information such as the number of triangles (NT), the number of vertices (NV), and
 pointers to each of the data dields. The metadata and subject data areas are not currently used, but the
 layout is designed to store embedded XML to allow a variety of information to be stored in the file.
 The triangle data is a vector of NT triangles. Each triangle consists of 3 32-bit integers, so this data block
 is NT * 3 * 4 bytes in length. Each triangle specifying a triple of vertex indices that point into the vertex 
 data block. The vertex indices are specified starting from zero with a maximum value of (NV-1). The vertex data
 block is of length NV * 3 * 4, as each vertex is represented by a triple of 32-bit floating point numbers.

 All fields after the triangle and vertex data are optional. Each data field has NV elements, with the size of
 the element depending on the particular attribute. For example, UV coordinates are specified as pairs of 32-bit
 floating point numbers, so the vertex UV block will be of length NV * 2 * 4 if it is present.

 The label and attribute fields can be used to store various types of scalar data. For example, the labels may be
 used to indicate anatomical areas or regions of interest. Similarly, the vertex attribute may be used to store
 curvature measures or activation measures.
*/


#include <vector>
#include <string>
#include <fstream>
#include <iostream>

namespace SILT
{
inline void endian_swap(unsigned short& x)
/// These endian swap functions were obtained from CodeGuru C++ FAQ:
/// http://www.codeguru.com/forum/showthread.php?t=292902
{
    x = (x>>8) | 
        (x<<8);
}

inline void endian_swap(unsigned int& x)
/// These endian swap functions were obtained from CodeGuru C++ FAQ:
/// http://www.codeguru.com/forum/showthread.php?t=292902
{
    x = (x>>24) | 
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
}

inline void endian_swap(int& x)
/// This function recasts an signed 16-bit integer to use the above unsigned integer byte-swapping
{
	endian_swap(*(unsigned int *)&x);
}


inline void endian_swap(float& x)
/// This function recasts a 32-bit floating point number to use the above unsigned integer byte-swapping
{
	endian_swap(*(unsigned int *)&x);
}


/**
 \brief  ByteOrder provides two simple functions for testing if a system is big endian or little endian.
 \author David Shattuck
*/
class ByteOrder
{
public:
	static bool bigEndian() ///< returns true if the current system stores data in big endian format
	{
		short int word = 0x0001;
		char *byte = (char *) &word;
		return (!byte[0]);
	}
	static bool littleEndian() ///< returns true if the current system stores data in little endian format
	{
		short int word = 0x0001;
		char *byte = (char *) &word;
		return (byte[0]!=0);
	}
};

/**
 \brief  BrainSuite Surface Format (dfs) file header
 \author David Shattuck

 The DFSHeader consists of:

 (1) a 12 byte header identifying the file type, data endianness, and header version number.<br>
 (2) a series of 4-byte integers specifying the size of the header, the numbers of triangles
     and vertices stored in the file, and the locations of other data.<br>

 The triangle and vertex arrays are required. All other fields are optional, and some of the
 fields specified in the header are now deprecated (e.g., stripSize and nStrips).  The field
 index represents the file position where the attribute vector can be read. If a field index
 is set to zero it indicates that the attributes for that field are not present in the file.
 */
class DFSHeader {
public:
	DFSHeader() : ///< The constructor for DFSHeader sets the headerType field and the headerSize. All other values are set to zero.
		headerSize(sizeof(DFSHeader)), metadataOffset(0), subjectDataOffset(0),
		nTriangles(0), nVertices(0), nStrips(0), stripSize(0),
		vertexNormalsOffset(0), vertexUVOffset(0), vertexColorsOffset(0),
		vertexLabelsOffset(0), vertexAttributesOffset(0)
	{
		headerType[0]='D';
		headerType[1]='F';
		headerType[2]='S';
		headerType[3]='_';
		headerType[4]=(ByteOrder::bigEndian()) ? 'B' : 'L';
		headerType[5]='E';
		headerType[6]=' ';
		headerType[7]='v';
		headerType[8]='2';
		headerType[9]='.';
		headerType[10]='0';
		headerType[11]='\0';
	}
	void swap() ///< performs byte-swapping on the fields in the header
	{
		endian_swap(headerSize);
		endian_swap(metadataOffset);
		endian_swap(subjectDataOffset);
		endian_swap(nTriangles);
		endian_swap(nVertices);
		endian_swap(nStrips);
		endian_swap(stripSize);
		endian_swap(vertexNormalsOffset);
		endian_swap(vertexUVOffset);
		endian_swap(vertexColorsOffset);
		endian_swap(vertexLabelsOffset);
		endian_swap(vertexAttributesOffset);
	}
	bool isSwapped() ///< returns true if the current contents of the header are in a different byte order than the system
	{
		if (headerType[4]=='B') { return ByteOrder::littleEndian(); }
		if (headerType[4]=='L') { return ByteOrder::bigEndian(); }
		return (headerSize>32000); // in case it is an old version of the header
	}
	char headerType[12];				///< 12 character string identifying the header type and format. This will be either DFS_BE v2.0 or DFS_LE v2.0 (null terminated in both cases)
	int headerSize;							///< size of complete header (i.e., offset of first data element)
	int metadataOffset;					///< start of metadata (not used presently)
	int subjectDataOffset;			///< start of subject data (not used presently)
	int nTriangles;							///< number of triangles
	int nVertices;							///< number of vertices
	int nStrips;								///< number of triangle strips (not used)
	int stripSize;							///< size of strip data (not used)
	int vertexNormalsOffset;		///< start of vertex normal data (0 if not in file)
	int vertexUVOffset;					///< start of surface parameterization data (0 if not in file)
	int vertexColorsOffset;			///< vertex color
	int vertexLabelsOffset;			///< vertex labels
	int vertexAttributesOffset;	///< vertex attributes (float32 array of length NV)
	char pad1[8*4*4-4];					///< this area previously contained a structure that is now deprecated; it shrinks when new fields are added to the header
};

template <class T>
inline std::istream &operator>>(std::istream &is, std::vector<T> &v) /// read a vector from an input stream (assumed binary)
{
	if (!v.empty())
		is.read((char *)&v[0],(std::streamsize)(v.size()*sizeof(T)));
	return is;
}

template <class T>
inline std::ostream &operator<<(std::ostream &os, std::vector<T> &v) /// read a vector to an output stream (assumed binary)
{
	if (!v.empty())
		os.write((char *)&v[0],(std::streamsize)(v.size()*sizeof(T)));
	return os;
}

/// \brief The Triangle class consists of 3 integer values that index into a vertex list.
/// \author David Shattuck
///
/// BrainSuite assumes that the list of vertices starts at zero.
class Triangle {
public:
	Triangle() : a(0), b(0), c(0) {}
	int a,b,c;
};


/// \brief The PointUV class consists of 2 32-bit floating point values that specify a 2D coordinate.
/// \author David Shattuck
///
/// The UV coordinates may be used for representing flattened surfaces or for performing texture
/// mapping onto a surface.
class PointUV {
public:
	PointUV() : u(0), v(0) {}
	float u,v;
};

/// \brief The Point3D class is specifies a 3D coordinate using 32-bit floating point values.
class Point3D {
public:
	Point3D() : x(0), y(0), z(0) {}
	float x,y,z;
};


/// \brief swap each element of a UV Point
inline void endian_swap(PointUV &p)
{
	endian_swap(p.u);
	endian_swap(p.v);
}


/// \brief swap each element of a triangle
inline void endian_swap(Triangle &t)
{
	endian_swap(t.a);
	endian_swap(t.b);
	endian_swap(t.c);
}


/// \brief swap each element of a 3D point
inline void endian_swap(Point3D &p)
{
	endian_swap(p.x);
	endian_swap(p.y);
	endian_swap(p.z);
}


/// \brief swap each element in the vector v
template <class T>
void swap_vec(std::vector<T> &v)
{
	const size_t n = v.size();
	for (size_t i=0;i<n;i++)
		endian_swap(v[i]);
}

/// \brief The DFSurface class contains several vectors of data, including triangle and vertex
/// arrays.
/// \author David Shattuck
///
/// The DFSurface class contains several vectors of data, including triangle and vertex
/// arrays. The DFSurface class also contains fields for optional data, including surface
/// normals, vertex colors, and other attributes.
class DFSurface {
public:	
	/// reads a file that is in DFS format; returns true if successful and false if unsuccessful.
	bool readDFS(std::string ifname) ///< the input filename (should end in .dfs).
	{
		std::ifstream ifile;
		ifile.open(ifname.c_str(), std::ios::binary);
		if (!ifile)
		{
			std::cerr<<"unable to open "<<ifname.c_str()<<std::endl;
			return false;
		}
		DFSHeader header;
		ifile.read((char *)&header,sizeof(header));
		bool swapped = header.isSwapped();
		if (swapped) header.swap();
		if (header.headerType[8]==1) // old headers had version number encoded in bytes 8-11
		{
			if (header.headerType[11]<2)
			{
				header.vertexLabelsOffset     = 0; // and versions before 1.0.0.2 didn't have labels
				header.vertexAttributesOffset = 0; // and versions before 1.0.0.2 didn't have attributes
			}
		}
		triangles.resize(header.nTriangles);
		vertices.resize(header.nVertices);
		vertexNormals.resize((header.vertexNormalsOffset>0)  ? header.nVertices : 0);
		vertexColors.resize((header.vertexColorsOffset>0) ? header.nVertices : 0);
		vertexUV.resize((header.vertexUVOffset>0) ? header.nVertices : 0);
		vertexLabels.resize((header.vertexLabelsOffset>0) ? header.nVertices : 0);
		vertexAttributes.resize((header.vertexAttributesOffset>0) ? header.nVertices : 0);

		ifile.seekg(header.headerSize+header.metadataOffset+header.subjectDataOffset,std::ios_base::beg);
		ifile>>triangles;
		ifile>>vertices;
		if (header.vertexNormalsOffset>0)
		{
			ifile.seekg(header.vertexNormalsOffset,std::ios_base::beg);
			ifile>>vertexNormals;
		}
		if (header.vertexColorsOffset>0)
		{
			ifile.seekg(header.vertexColorsOffset,std::ios_base::beg);
			ifile>>vertexColors;
		}
		if (header.vertexUVOffset>0)
		{
			ifile.seekg(header.vertexUVOffset,std::ios_base::beg);
			ifile>>vertexUV;
		}
		if (header.vertexLabelsOffset>0)
		{
			ifile.seekg(header.vertexLabelsOffset,std::ios_base::beg);
			ifile>>vertexLabels;
		}
		if (header.vertexAttributesOffset>0)
		{
			ifile.seekg(header.vertexAttributesOffset,std::ios_base::beg);
			ifile>>vertexAttributes;
		}
		if (swapped)
		{
			swap_vec(triangles);
			swap_vec(vertices);
			swap_vec(vertexNormals);
			swap_vec(vertexColors);
			swap_vec(vertexUV);
			swap_vec(vertexLabels);
			swap_vec(vertexAttributes);
		}
		return true;
	}
/// writes a surface object as a DFS file; returns true if successful and false if unsuccessful.
	bool writeDFS(std::string ofname)	///< the output filename (should end in .dfs).
	{
		std::ofstream ofile;
		ofile.open(ofname.c_str(), std::ios::binary);
		if (!ofile) return false;
		DFSHeader header;
		const int nv = (int)vertices.size();
		const int nt = (int)triangles.size();
		header.nTriangles = nt;
		header.nVertices = nv;
		int pos = sizeof(header) + sizeof(Triangle)*nt + sizeof(Point3D)*nv;
		if (vertexNormals.size()==nv)    { header.vertexNormalsOffset = pos; pos += sizeof(vertexNormals[0])*nv; }
		if (vertexColors.size()==nv)     { header.vertexColorsOffset = pos; pos += sizeof(vertexColors[0])*nv; }
		if (vertexUV.size()==nv)		     { header.vertexUVOffset = pos; pos += sizeof(vertexUV[0])*nv; }
		if (vertexLabels.size()==nv)		 { header.vertexLabelsOffset = pos; pos += sizeof(vertexLabels[0])*nv; }
		if (vertexAttributes.size()==nv) { header.vertexAttributesOffset = pos; pos += sizeof(vertexAttributes[0])*nv; }

		ofile.write((char *)&header,sizeof(header));
		ofile<<triangles;
		ofile<<vertices;
		if (vertexNormals.size()==nv) ofile<<vertexNormals;
		if (vertexColors.size()==nv) ofile<<vertexColors;
		if (vertexUV.size()==nv) ofile<<vertexUV;
		if (vertexLabels.size()==nv) ofile<<vertexLabels;
		if (vertexAttributes.size()==nv) ofile<<vertexAttributes;
		return true;
	}
	std::vector<Triangle> triangles; ///< triangle index array; indices start at zero
	std::vector<Point3D> vertices;   ///< array of mesh vertices
	std::vector<Point3D> vertexNormals; ///< array of vertex normals (optional)
	std::vector<Point3D> vertexColors; ///< array of vertex colors (optional)
	std::vector<PointUV> vertexUV; ///< array of surface-based coordinates for flat-mapping (optional)
	std::vector<unsigned short> vertexLabels; ///< array of labels for the vertices (optional)
	std::vector<float> vertexAttributes; ///< array of attributes for the vertices (optional)
};

} // end of namespace SILT

#endif

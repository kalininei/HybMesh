#ifndef HMCONT2D_ALGOS_HPP
#define HMCONT2D_ALGOS_HPP

namespace HMCont2D{

//Offset contour
enum class OffsetTp{
	CLOSED_POLY,
	OPEN_ROUND,
	OPEN_BUTT,
};

enum class PartitionTp{
	IGNORE_ALL,
	KEEP_ALL,
	KEEP_SHAPE
};

//structs declaration to make contour/tree see each other
struct Contour;
struct ContourTree;
//

};


#endif

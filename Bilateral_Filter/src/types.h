#ifndef TYPES_H
#define TYPES_H

#include<list>
#include<vector>
#include<set>
#include<map>
#include<deque>
#include "Point.h"

class Sample;


typedef std::deque<Sample> Sample_deque;
typedef std::deque<Sample*> Sample_star_deque;

typedef std::map<double, Sample*> Neighbor_star_map;
typedef std::list<Sample*> Neighbor_star_list;

typedef std::list<double> Distance_list;
typedef std::list<Sample*> Sample_star_list;

#include "Octree.h"
typedef TOctree<Sample> Octree;

#include "OctreeNode.h"
typedef TOctreeNode<Sample> OctreeNode;

#include "OctreeIterator.h"
typedef TOctreeIterator<Sample> OctreeIterator;


#include "BilateralFilter.h"
typedef TBilateralFilter<Sample> BilateralFilter;


#endif

/**
 * @brief defines an octree node
 */

#ifndef OCTREENODE_H
#define OCTREENODE_H

#include <cstdlib>
#include <deque>
#include <iostream>
#include <fstream>
#include <cassert>

#include "Point.h"
/**
 * @class TOctreeNode
 * @brief Implements a generic node for a generic octree
 * Templated class implementing a node of the octree, leaf nodes
 * contain the input points.
 */
template<class T>
class TOctreeNode
{
    protected :
        TOctreeNode<T> *m_parent;

        /** @brief pointer to the eight children of the node
         */
        TOctreeNode<T> *m_child[8];

        /** @brief 
         * child number of the node (depends on the relative location
         * of the node to the middle of its parent)
         * \verbatim	
         0 *-------4
         /|      /|
         2-------6 |
         | 1-----|-5
         |/      |/
         3-------7
         axis:
         x: along direction 0->4
         y: along direction 0->2
         z: along direction 0->1
         \endverbatim
         * REMARK: this convention is not important because we are dealing with 
         * a cube.
         * It will not affect neither representation nor computations
         */
        unsigned int m_nchild;

        /** @brief number of points included in the node or in the
         * node's children 
         */
        unsigned int m_npts;

        /** @brief origin of the node*/
        Point m_origin;

        /** @brief level of the node*/
        unsigned int m_depth;

        /**
         * @brief x locational code
         */
        unsigned int m_xloc;

        /**
         * @brief y locational code
         */
        unsigned int m_yloc;

        /**
         * @brief z locational code
         */
        unsigned int m_zloc;

        /** @brief size of the node side*/
        double m_size;

        /** @brief list of points contained in the node
         * (empty if the node is not a leaf)
         * the first list are the input points contained in this node
         * the second and third lists are the filtered sets
         * (2 lists are necessary for buffering)
         */
        std::deque<T> m_points[3];

    public :
        /**
         *  @brief Default constructor initializes all variables
         */
        TOctreeNode();
        TOctreeNode(Point & origin, double size, unsigned int depth);
        ~TOctreeNode();

        /**
         * @brief set node size
         * @param size desired node size
         */
        void setSize (double size);

        /**
         * @brief get node size
         * @return size of the node
         */
        double getSize() const;

        /**@brief returns the number of points contained in the node
         * @param index index of the set to process
         * @return number of points contained in the node (if leaf), 0 otherwise
         */
        unsigned int getNpts(unsigned int index) const;

        /**
         * @brief set child number (depends on the relative location of
         * the node to the middle of its parent)
         * @param a child number
         */
        void setNchild(unsigned int a);

        /** @brief get the child number
         * @return child number 
         */
        const unsigned int& getNChild() const;

        /**
         * @brief set the origin of the node by a point structure
         * @param pt Origin point
         */
        void setOrigin (Point &pt);

        /** @brief get the origin of a given node
         * @return origin
         */
        Point getOrigin() const;

        /**
         * @brief set parent of a node
         * @param node pointer to the parent node
         */
        void setParent ( TOctreeNode<T> *node );
        TOctreeNode<T>* getParent() const;

        /** @brief get child of a node
         * @param index of the child
         */
        TOctreeNode<T>* getChild(unsigned int index);
        /**
         * @brief set level of a node
         */
        void setDepth (unsigned int l );
        unsigned int getDepth() const;
        /**
         * @brief check if a point given by its coordinates is inside a node
         */
        bool isInside(double x,double y,double z) const;
        bool isInside(const Point &p) const;
        /**
         * @brief check if a point given by its coordinates is inside or
         * in a band around the node
         */
        bool isInside(const Point &p, double d) const;

        unsigned int getXLoc() const;
        unsigned int getYLoc() const;
        unsigned int getZLoc() const;

        /**
         * @brief locationnal code method
         */
        void setXLoc(unsigned int Xloc);
        void setYLoc(unsigned int Yloc);
        void setZLoc(unsigned int Zloc);

        /** @brief get a pointer to the list of points
         * @param index index of the point list
         * @return pointer to the beginning of the list
         */
        typename std::deque<T>::iterator points_begin(unsigned int index);

        /** @brief get a pointer to the end of the list of points
         * @param index index of the point list
         * @return pointer to the end of 'points'
         */
        typename std::deque<T>::iterator points_end(unsigned int index);

        /** @brief get a const pointer to the list of points
         * @param index index of the point list
         * @return const pointer to the beginning of the list
         */
        typename std::deque<T>::const_iterator points_begin(
                                              unsigned int index) const;

        /** @brief get a const pointer to the end of the list of points
         * @param index index of the point list
         * @return const pointer to the end of 'points'
         */
        typename std::deque<T>::const_iterator points_end(
                                               unsigned int index) const;

        /** @brief add a point to the list of points included in the cell
         * PREREQUISITE: the node is a leaf in the octree
         * @param pt point to add
         */
        void addInitialPoint(T &pt);

        /** @brief add a point to the list of points included in the cell
         * PREREQUISITE: the node is a leaf in the octree
         * @param pt point to add
         * @param index index of the point list
         */
        void addPoint(T &pt,unsigned int index);

        /** @brief build the i^th child of the node
         * @param index child index
         * @param origin origin of the node
         * @return pointer to the created node 
         */
        TOctreeNode<T>* initializeChild(unsigned int index, Point origin);

        /** @brief clear point sets in all the children of the node
         * (and itself) 
         * @param index index of the set to clear
         */
        void clearSet(unsigned int index);
};

template<class T>
TOctreeNode<T>::TOctreeNode()
{
    for(int i = 0 ; i <8 ; i++)
        m_child[i] = NULL;
    m_parent = NULL;
    m_xloc = m_yloc = m_zloc =0;
    m_depth = 0;
    m_npts = 0;
    m_origin = Point();
    m_size = 0.0;
}

template<class T>
TOctreeNode<T>::TOctreeNode(Point& origin, double size, unsigned int depth)
{
    for(int i = 0 ; i <8 ; i++)
        m_child[i] = NULL;
    m_parent = NULL;
    m_xloc = m_yloc = m_zloc =0;
    m_depth = depth;
    m_npts = 0;
    m_origin = origin;
    m_size = size;
}

template<class T>
TOctreeNode<T>::~TOctreeNode()
{
    for(int i = 0; i < 3 ; ++i)
        m_points[i].clear();
    
    m_xloc = m_yloc = m_zloc =0;
    m_depth = 0;
    m_npts = 0;
    for(int i = 0; i<8 ; i++)
        delete m_child[i];
    m_parent = NULL;
    m_origin = Point();
    m_size = 0.0;
}


template<class T>
unsigned int TOctreeNode<T>::getDepth() const
{
    return m_depth;
}

template<class T>
void TOctreeNode<T>::setDepth(unsigned int l)
{
    m_depth = l;
}

template<class T>
double TOctreeNode<T>::getSize() const
{
    return m_size;
}

template<class T>
void TOctreeNode<T>::setSize(double size)
{
    m_size = size;
}


template<class T>
unsigned int TOctreeNode<T>::getNpts(unsigned int index) const
{
    return m_points[index].size();
}


template<class T>
const unsigned int& TOctreeNode<T>::getNChild() const
{
    return m_nchild;
}

template<class T>
void TOctreeNode<T>::setNchild(unsigned int a)
{
    m_nchild = a;
}


template<class T>
void TOctreeNode<T>::setParent(TOctreeNode *parent)
{
    m_parent = parent;
}

template<class T>
TOctreeNode<T>* TOctreeNode<T>::getParent() const
{
    return m_parent;
}

template<class T>
TOctreeNode<T>* TOctreeNode<T>::getChild(unsigned int index) 
{
    unsigned int i = index % 8;
    return m_child[i];
}


template<class T>
unsigned int TOctreeNode<T>::getXLoc() const
{
    return m_xloc;
}

template<class T>
void TOctreeNode<T>::setXLoc(unsigned int xloc)
{
    m_xloc = xloc;
}

template<class T>
unsigned int TOctreeNode<T>::getYLoc() const
{
    return m_yloc;
}

template<class T>
void TOctreeNode<T>::setYLoc(unsigned int yloc)
{
    m_yloc = yloc;
}

template<class T>
unsigned int TOctreeNode<T>::getZLoc() const
{
    return m_zloc;
}

template<class T>
void TOctreeNode<T>::setZLoc(unsigned int zloc)
{
    m_zloc = zloc;
}

template<class T>
bool TOctreeNode<T>::isInside(double x, double y, double z) const
{
    
    if ( ( x >= m_origin.x() )&&( x < m_origin.x() + m_size)
        && ( y >= m_origin.y() )&&( y < m_origin.y() + m_size)
        && ( z >= m_origin.z() )&&( z < m_origin.z() + m_size))
        return true;
    else
        return false;
}

template<class T>
bool TOctreeNode<T>::isInside(const Point &p) const
{
    if ( ( p.x() >= m_origin.x() )&&( p.x() < m_origin.x() + m_size)
        && ( p.y() >= m_origin.y() )&&( p.y() < m_origin.y() + m_size)
        && ( p.z() >= m_origin.z() )&&( p.z() < m_origin.z() + m_size))
        return true;
    else
        return false;
}

template<class T>
bool TOctreeNode<T>::isInside(const Point &p, double d) const
{
    double offset = m_size + d;
    if ( ( p.x() >= m_origin.x() - d )&&( p.x() < m_origin.x() + offset)
        && ( p.y() >= m_origin.y() - d )&&( p.y() < m_origin.y() + offset)
        && ( p.z() >= m_origin.z() - d )&&( p.z() < m_origin.z() + offset))
        return true;
    else
        return false;
}

template<class T>
void TOctreeNode<T>::setOrigin(Point& pt)
{
    m_origin = pt;
}


template<class T>
Point TOctreeNode<T>::getOrigin() const
{
    return m_origin;
}

template<class T>
typename std::deque<T>::iterator TOctreeNode<T>::points_begin(unsigned int 
index)
{
    return m_points[index].begin();
}

template<class T>
typename std::deque<T>::iterator TOctreeNode<T>::points_end(unsigned int index)
{
    return m_points[index].end();
}

template<class T>
typename std::deque<T>::const_iterator TOctreeNode<T>::points_begin(
    unsigned int index) const
{
    return m_points[index].begin();
}

template<class T>
typename std::deque<T>::const_iterator TOctreeNode<T>::points_end(
unsigned int index) const
{
    return m_points[index].end();
}

template<class T>
void TOctreeNode<T>::addInitialPoint(T &t)
{
    m_points[0].push_back(t);
    m_npts++;
}

template<class T>
void TOctreeNode<T>::addPoint(T &t, unsigned int index)
{
    m_points[index].push_back(t);
}

template<class T>
TOctreeNode< T >* TOctreeNode<T>::initializeChild(unsigned int index, Point origin){
    double size = m_size/2.0;
    unsigned int depth= m_depth -1;
    m_child[index] = new TOctreeNode<T>(origin, size, depth);
    m_child[index]->setParent(this);
    m_child[index]->setNchild(index);

    return m_child[index];
}

template<class T>
void TOctreeNode<T>::clearSet(unsigned int index){
    if(getDepth() == 0){
        m_points[index].clear();
    }
    else{
        for(unsigned int i=0; i<8; ++i)
            if(getChild(i) != NULL)
                getChild(i)->clearSet(index);
    }
}
#endif
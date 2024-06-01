
#ifndef FILEIO_H
#define FILEIO_H

#include <fstream>
#include <iostream>

#include "Octree.h"
#include "Sample.h"
#include "types.h"


/**
 * @class FileIO
 * @brief class providing access to input/output operations
 *
 * This class defines methods for reading points from an input file
 * and saving the resulting mesh in a ply file.
 */
class FileIO
{
   public :

   /** @brief constructor*/
   FileIO();

   /** @brief destructor*/
   ~FileIO();

   public :

   /** @brief estimate the radius
    * heuristic to find the radius (see paper for details)
    * @param filename name of the point set to estimate the radius from
    * @return a radius estimation
    */
   static double estimateRadius(const char *filename);

   /** @brief read points from a file
    * @param filename name of the file to read points from
    * @param octree to sort and store the points in
    * @param min_radius if positive, create the octree such that the smallest
    * cell has size 2*min_radius
    * @return true if saving was successful
    */
   static bool readAndSortUnorientedPoints(const char *filename, Octree &octree, 
                                 double min_radius = -1);

   /** @brief read points from a file
    * @param filename name of the file to read points from
    * @param octree to sort and store the points in
    * @param min_radius if positive, create the octree such that the smallest
    * cell has size 2*min_radius
    * @return true if saving was successful
    */
   static bool readAndSortOrientedPoints(const char *filename, Octree &octree, 
                                 double min_radius = -1);
   /** @brief save points from an octree to a file
    * @param filename name of the file to save to
    * @param octree octree to save the points from
    * @param index index of the octree set to save from
    * @param isoriented should be set to true if the input point set was
    * oriented, false otherwise.
    * @return true if saving was successful
    */
   static bool savePoints(const char* filename, Octree &octree,
                          unsigned int index, bool isoriented);

    private :

    /** @brief save all vertices contained in a node
     * @param node node to save from
     * @param index index of the octree set to save from
     * @param f stream to save to
    * @param isoriented should be set to true if the input point set was
    * oriented, false otherwise.
     */
    static void saveContent(OctreeNode *node, Octree &octree, unsigned int index,
                            std::ofstream &f, bool isoriented);
};


#endif

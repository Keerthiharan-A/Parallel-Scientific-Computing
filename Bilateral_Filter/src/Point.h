/**
 * @file Point.h
 * @brief Declares a generic point class
 */ 

#ifndef POINT_H
#define POINT_H

//The most standard 3D point structure: only contains three coordinates

class Point
{
    public :
        Point();

        /** @param x x coordinate
         * @param y y coordinate
         * @param z z coordinate
         **/
        Point(double x, double y, double z);


        /** @brief access x coordinate
         * @return x
         */
        double x() const;
        
        double y() const;

        double z() const;


    private :
        
        /** @brief x coordinate*/
        double m_x;
        
        double m_y;
        
        double m_z;

};

#endif
/** @brief Declaration of a vertex object 
 */

#ifndef SAMPLE_H
#define SAMPLE_H

#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <set>

#include "Point.h"
#include "types.h"


/**
 * @class Sample
 * @brief Input samples to be filtered
 * 
 */
class Sample : public Point
{
    /** @brief overloading operator <<
     * @param out output stream
     * @param v vertex
     * @return output stream
     */
    friend std::ostream& operator << (std::ostream& out, const Sample& v);
  
    private : //properties
    
        /** @brief nx, ny, nz normal coordinates*/
        double m_nx, m_ny, m_nz;
        
        /** @brief sample index */
        int m_index;
        

    public : //constructor+destructor
    
        /** @brief default constructor*/
        Sample();
        
        /** @brief constructor from coordinates and normal*/
        Sample(double x, double y, double z);
        
        /** @brief constructor from coordinates and normal*/
        Sample(double x, double y, double z, double nx, double ny, double nz);
          
        /** @brief default destrictor*/
        ~Sample();
    
    public :

       //accessors + modifiers
        /** @brief get index of the vertex
         @return index
         */
        int index() const;
        
        /** @brief set the vertex index
         * @param index to set
         */
        void setIndex(int index);

        /** @brief get x normal component
         * @return x normal component nx
         */
        double nx() const;
        
        /** @brief get y normal component
         * @return y normal component ny
         */
        double ny() const;
        
        /** @brief get z normal component
         * @return z normal component nz
         * */
        double nz() const;
};



#endif

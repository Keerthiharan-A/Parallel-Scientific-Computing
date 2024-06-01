/**
 * @brief defines generic class Point
 */
#include "Point.h"
#include <cstdlib>

Point::Point()
{
    m_x = m_y = m_z = 0;
}
Point::Point(double x, double y, double z)
{
    m_x = x;
    m_y = y;
    m_z = z;
}

double Point::x() const
{
    return m_x;
}
double Point::y() const
{
    return m_y;
}
double Point::z() const
{
    return m_z;
}
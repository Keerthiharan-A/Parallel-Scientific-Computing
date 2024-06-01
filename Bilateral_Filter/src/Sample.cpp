#include "Sample.h"
#include "utilities.h"
#include <cstdio>
#include <iostream>

Sample::Sample() : Point()
{
    m_nx=m_ny=m_nz=0.0;
    m_index = -1;
}

Sample::Sample(double x, double y, double z, double nx, double ny, double nz)
                : Point(x,y,z)
{
    m_nx = nx;
    m_ny = ny;
    m_nz = nz;
    m_index = -1;
}

Sample::Sample(double x, double y, double z)
                : Point(x,y,z)
{
    m_nx = 0;
    m_ny = 0;
    m_nz = 0;
    m_index = -1;
}

Sample::~Sample()
{
    m_nx=m_ny=m_nz=0.0;
    m_index = -1;
}

int Sample::index() const
{
    return m_index;
}

void Sample::setIndex(int index)
{
    m_index = index;
}


double Sample::nx() const
{
    return m_nx;
}

double Sample::ny() const
{
    return m_ny;
}

double Sample::nz() const
{
    return m_nz;
}



std::ostream& operator << (std::ostream& out, const Sample& v)
{
    out << v.x() << "\t" << v.y() << "\t" << v.z()
        << "\t" << v.nx() << "\t" << v.ny() << "\t" << v.nz();
    return out;
}




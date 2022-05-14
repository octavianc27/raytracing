#include <cmath>

#include "Vector3f.h"

Vector3f::Vector3f() : Vector3f(0, 0, 0)
{
}

Vector3f::Vector3f(const float x, const float y, const float z)
{
    setValues(x, y, z);
}

float
Vector3f::getX() const
{
    return x;
}

float
Vector3f::getY() const
{
    return y;
}

float
Vector3f::getZ() const
{
    return z;
}

void
Vector3f::setX(const float value)
{
    x = value;
}

void
Vector3f::setY(const float value)
{
    y = value;
}

void
Vector3f::setZ(const float value)
{
    z = value;
}

void
Vector3f::setValues(const float x, const float y, const float z)
{
    setX(x);
    setY(y);
    setZ(z);
}

Vector3f
Vector3f::operator + (const Vector3f & vector) const
{
    return Vector3f(x + vector.getX(), y + vector.getY(), z + vector.getZ());
}

Vector3f
Vector3f::operator - (const Vector3f & vector) const
{
    return Vector3f(x - vector.getX(), y - vector.getY(), z - vector.getZ());
}

Vector3f
Vector3f::operator * (const float scale) const
{
    return Vector3f(x * scale, y * scale, z * scale);
}

Vector3f
Vector3f::operator / (const float scale) const
{
    return Vector3f(x / scale, y / scale, z / scale);
}

Vector3f
Vector3f::crossProduct(const Vector3f & vector) const
{
    return Vector3f(y * vector.getZ() - z * vector.getY(),
                    z * vector.getX() - x * vector.getZ(),
                    x * vector.getY() - y * vector.getX());
}

float
Vector3f::innerProduct(const Vector3f & vector) const
{
    return (x * vector.getX() + y * vector.getY() + z * vector.getZ());
}

float
Vector3f::length() const
{
    return sqrt(x*x + y*y + z*z);
}

float
Vector3f::length_squared() const
{
    return (x*x + y*y + z*z);
}

Vector3f
Vector3f::operator - () const
{
    return Vector3f(-x, -y, -z);
}

void
Vector3f::opposite()
{
    setValues(-x, -y, -z);
}

Vector3f
Vector3f::normalized() const
{
    const float distance = length();

    return Vector3f(x / distance, y / distance, z / distance);
}

void
Vector3f::normalize()
{
    const float distance = length();

    setValues(x / distance, y / distance, z / distance);
}
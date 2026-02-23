#include "polygon2D.h"
#include <cmath>
#include <vector>
#include <algorithm>

using namespace std;

double polygon2D::polyArea2d(vector<pair<double, double>> & vertices)
{
    vector<pair<double,double>> sortedVertices = sortVertices(vertices);
    double area = 0.0;
    for (int i = 0; i < sortedVertices.size() ; i++)
    {
        area += sortedVertices[i].first * sortedVertices[(i + 1) % sortedVertices.size()].second - sortedVertices[(i + 1) % sortedVertices.size()].first * sortedVertices[i].second;
    }
    return 0.5 * abs(area);
    
}

double polygon2D::polyPerimeter2d(vector<pair<double, double>> & vertices)
{
    vector<pair<double,double>> sortedVertices = sortVertices(vertices);
    double perimeter = 0.0;
    for (int i = 0; i < sortedVertices.size(); i++)
    {
        perimeter += sqrt(pow((sortedVertices[i].first - sortedVertices[(i+1)%sortedVertices.size()].first),2) + pow((sortedVertices[i].second - sortedVertices[(i+1)%sortedVertices.size()].second),2));
    }
    
    return perimeter;
}

vector<double> polygon2D::polyAngle2d(vector<pair<double,double>> & vertices)
{
    vector<pair<double,double>> sortedVertices = sortVertices(vertices);
    vector<double> sortedAngle;
    vector<double> outputAngle;
    for(int i=0;i<sortedVertices.size();i++)
    {
        double ang;
        int preI;
        int aftI;
        if(i==0)
        {
            preI = sortedVertices.size()-1;
            aftI = i + 1;
        }
        else if(i==sortedVertices.size()-1)
        {
            preI = i - 1;
            aftI = 0;
        }
        else
        {
            preI = i-1;
            aftI = i+1;
        }

        ang = calculateAngle(sortedVertices[preI],sortedVertices[i],sortedVertices[aftI]);
        sortedAngle.push_back(ang);
    }


    for(int i=0;i<vertices.size();i++)
    {
        for(int j = 0; j < sortedVertices.size();j++)
        {
            if(vertices[i] == sortedVertices[j])
            {
                outputAngle.push_back(sortedAngle[j]);
            }
        }
    }

    return outputAngle;
}


pair<double, double> polygon2D::getCentroid2d(const vector<pair<double, double>> & vertices)
{
    double cx = 0.0;
    double cy = 0.0;
    for (const auto &vertex : vertices)
    {
        cx += vertex.first;
        cy += vertex.second;
    }
    cx /= vertices.size();
    cy /= vertices.size();

    return {cx, cy};
}

pair<double, double> polygon2D::getCentroid2d(const vector<tuple<double, double,int>> & vertices)
{
    double cx = 0.0;
    double cy = 0.0;
    for (const auto &vertex : vertices)
    {
        cx += get<0>(vertex);
        cy += get<1>(vertex);
    }
    cx /= vertices.size();
    cy /= vertices.size();

    return {cx, cy};
}



double polygon2D::calculateAngle(const pair<double, double> & vertex1, const pair<double, double> & vertex2, const pair<double, double> & vertex3)
{
    pair<double,double> AB = {{vertex1.first-vertex2.first},{vertex1.second - vertex2.second}};
    pair<double,double> AC = {{vertex3.first-vertex2.first},{vertex3.second - vertex2.second}};
    double dot_prod = AB.first * AC.first + AB.second * AC.second;
    double mag_AB = sqrt(AB.first * AB.first + AB.second * AB.second);
    double mag_AC = sqrt(AC.first * AC.first + AC.second * AC.second);

    double angle = acos(dot_prod/(mag_AB*mag_AC)) * 180/ M_PI;
    return angle;
}

bool polygon2D::compareAngles(const pair<double, double> & vertex1, const pair<double, double> & vertex2, const pair<double, double> & centroid)
{
    double angle1 = atan2(vertex1.second - centroid.second, vertex1.first - centroid.first);
    double angle2 = atan2(vertex2.second - centroid.second, vertex2.first - centroid.first);
    return angle1 < angle2;
}

bool polygon2D::compareAngles(const tuple<double, double,int> & vertex1, const tuple<double, double,int> & vertex2, const pair<double, double> & centroid)
{
    double angle1 = atan2(get<1>(vertex1) - centroid.second, get<0>(vertex1) - centroid.first);
    double angle2 = atan2(get<1>(vertex2) - centroid.second, get<0>(vertex2) - centroid.first);
    return angle1 < angle2;
}

vector<pair<double, double>> polygon2D::sortVertices(vector<pair<double, double>> & vertices)
{
    pair<double,double> centroid = getCentroid2d(vertices);
    sort(vertices.begin(), vertices.end(), [&](const pair<double, double> &a, const pair<double, double> &b)
       { return compareAngles(a, b, centroid); });
    return vertices;
}

template <typename T>
vector<tuple<double, double, T>> polygon2D::sortVertices(vector<tuple<double, double, T>> &vertices)
{
  pair<double, double> centroid = getCentroid2d(vertices);
  sort(vertices.begin(), vertices.end(),
       [&](const tuple<double, double, T> &a, const tuple<double, double, T> &b) {
         return compareAngles(a, b, centroid);
       });
  return vertices;
}

vector<tuple<double, double,int>> polygon2D::sortVertices(vector<tuple<double, double,int>> & vertices)
{
    pair<double,double> centroid = getCentroid2d(vertices);
    sort(vertices.begin(), vertices.end(), [&](const tuple<double, double,int> &a, const tuple<double, double,int> &b)
       { return compareAngles(a, b, centroid); });
    return vertices;
}

vector<tuple<double, double,int>> polygon2D::sortVerticesPartly(vector<tuple<double, double,int>> & vertices, const pair<double, double> & centroid)
{
    sort(vertices.begin(), vertices.end(), [&](const tuple<double, double,int> &a, const tuple<double, double,int> &b)
       { return compareAngles(a, b, centroid); });
    return vertices;
}

template <typename T>
vector<tuple<double, double, T>> polygon2D::sortVerticesPartly(vector<tuple<double, double, T>> & vertices, const pair<double, double> & centroid)
{
    sort(vertices.begin(), vertices.end(), [&](const tuple<double, double, T> &a, const tuple<double, double, T> &b)
       { return compareAngles(a, b, centroid); });
    return vertices;
}

bool polygon2D::isPointInsidePolygon(const pair<double,double> &point, const vector<pair<double,double>> &vertices) 
{
    int count = 0;

    double minX = vertices[0].first;
    double maxX = vertices[0].first;
    double minY = vertices[0].second;

    
    // Find the minimum and maximum x and y coordinates of the polygon
    for (int i = 1; i < vertices.size(); i++) {
        minX = min(minX, vertices[i].first);
        maxX = max(maxX, vertices[i].first);
        minY = min(minY, vertices[i].second);
    }
    
    // Check if the point is outside the bounding box of the polygon
    if (point.first < minX || point.first > maxX || point.second < minY) {
        return false;
    }
    
    // Loop through each edge of the polygon
    for (int i = 0; i < vertices.size(); i++) {
        int j = (i + 1) % vertices.size(); // Next vertex
        
        // Check if the line segment from 'polygon[i]' to 'polygon[j]' intersects with the ray from 'p' horizontally to the right
        if ((vertices[i].second > point.second) != (vertices[j].second > point.second)) { // The point is between the y-coordinates of the edge
            // Compute the x-coordinate of the intersection
            double xIntersect = (vertices[j].first - vertices[i].first) * (point.second - vertices[i].second) / (vertices[j].second - vertices[i].second) + vertices[i].first;
            if (point.first < xIntersect) { // The intersection point is to the right of 'p'
                count++;
            }
        }
    }
    
    // 'count' is odd if the point is inside the polygon
    return count % 2 == 1;
}

vector<double> polygon2D::computeDistanceToEdgeCenter(const std::pair<double, double> & point,
                                                      const std::vector<pair<double, double>> & sortedVertices)
{
    vector<double> dist;

    for(int i = 0; i < sortedVertices.size();++i)
    {
        int j = (i + 1) % sortedVertices.size();
        
        double cx = 0.5 * (sortedVertices[i].first + sortedVertices[j].first);
        double cy = 0.5 * (sortedVertices[i].second + sortedVertices[j].second);
        dist.push_back(sqrt(pow(point.first-cx,2) + pow(point.second - cy,2)));
    }

  return dist;
}

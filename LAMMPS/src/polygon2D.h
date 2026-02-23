#ifndef POLYGON2D_H 
#define POLYGON2D_H
#include <vector>
#include <tuple>

class polygon2D
{
public:
    static double polyArea2d(std::vector<std::pair<double, double>> &);
    static double polyPerimeter2d(std::vector<std::pair<double, double>> &);
    static std::vector<double> polyAngle2d(std::vector<std::pair<double, double>> &);
    static std::pair<double, double> getCentroid2d(const std::vector<std::pair<double, double>> &);
    static std::pair<double, double> getCentroid2d(const std::vector<std::tuple<double, double,int>> &);
    static double calculateAngle(const std::pair<double, double> &,const std::pair<double, double> &,const std::pair<double, double> &);

    static std::vector<std::pair<double, double>> sortVertices(std::vector<std::pair<double, double>> &);
    static std::vector<std::tuple<double, double,int>> sortVertices(std::vector<std::tuple<double, double,int>> &);
    
    template <typename T>
    static std::vector<std::tuple<double, double, T>> sortVertices(std::vector<std::tuple<double, double, T>> &);
    static std::vector<std::tuple<double, double,int>> sortVerticesPartly(std::vector<std::tuple<double, double,int>> &, const std::pair<double, double> &);
    template <typename T>
    static std::vector<std::tuple<double, double, T>> sortVerticesPartly(std::vector<std::tuple<double, double, T>> &, const std::pair<double, double> &);
    static bool compareAngles(const std::pair<double, double> &, const std::pair<double, double> &, const std::pair<double, double> &);
    
    static bool compareAngles(const std::tuple<double, double,int> &, const std::tuple<double, double,int> &, const std::pair<double, double> &);

    static bool isPointInsidePolygon(const std::pair<double,double>&, const std::vector<std::pair<double,double>> &);

    static std::pair<double,double> computeCircumCenter(const std::pair<double,double>&, const std::pair<double,double>&,const std::pair<double,double>&);
    static std::vector<double> computeDistanceToEdgeCenter(const std::pair<double,double>&, const std::vector<std::pair<double,double>> &);
};

#endif
#include <iostream>
#include <math.h>
#include <fstream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <stdexcept>
#include <sstream>
#include <Eigen/Dense>
#include <iomanip>

class coordinate
{
public:
    double x;
    double y;
    double z;
    coordinate();
    coordinate(double x, double y, double z);
    ~coordinate();

coordinate operator + (coordinate const &obj)
{
    coordinate res;
    res.x = this->x + obj.x;
    res.y = this->y + obj.y;
    res.z = this->z + obj.z;
    return res;
}

coordinate operator - (coordinate const &obj)
{
    coordinate res;
    res.x = this->x - obj.x;
    res.y = this->y - obj.y;
    res.z = this->z - obj.z;
    return res;
}

coordinate operator / (coordinate const &obj)
{
    coordinate res;
    res.x = this->x / obj.x;
    res.y = this->y / obj.y;
    res.z = this->z / obj.z;
    return res;
}

void operator = (coordinate const &obj)
{
    x = obj.x;
    y = obj.y;
    z = obj.z;

}
};

coordinate::coordinate()
{
    this-> x = 0;
    this-> y = 0;
    this-> z = 0;
}

coordinate::coordinate(double x, double y, double z)
{
    this-> x = x;
    this-> y = y;
    this-> z = z;
}

coordinate::~coordinate()
{
}

double distance (coordinate a, coordinate b)
{
    return std::sqrt((a.x - b.x)*(a.x - b.x) + (a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
}


class ray
{
public:
    int tag;
    std::vector <double> coeff_x;
    std::vector <double> coeff_y;
    std::vector <double> coeff_z;
    ray();
    ray(coordinate a, coordinate b);
    void SetCoeff(coordinate a, coordinate b)
    {
    this-> coeff_x[0] = a.x;
    this-> coeff_y[0] = a.y;
    this-> coeff_z[0] = a.z;
    this-> coeff_x[1] = -b.x+a.x;
    this-> coeff_y[1] = -b.y+a.y;
    this-> coeff_z[1] = -b.z+a.z;
    }
    ~ray();
    
};

ray::ray()
{
    std::vector <double> coeff_x (2,0);
    std::vector <double> coeff_y (2,0);
    std::vector <double> coeff_z (2,0);
    this-> tag = 0;
    this-> coeff_x = coeff_x;
    this-> coeff_y = coeff_y;
    this-> coeff_z = coeff_z;
}

ray::ray(coordinate a, coordinate b)
{
    std::vector <double> coeff_x = {a.x,-b.x+a.x};
    std::vector <double> coeff_y = {a.y,-b.y+a.y};
    std::vector <double> coeff_z = {a.z,-b.z+a.z};
    //std::cout << "ray = " << a.z << a.z-b.z << std::endl;
    this-> tag = 0;
    this-> coeff_x = coeff_x;
    this-> coeff_y = coeff_y;
    this-> coeff_z = coeff_z;
}

ray::~ray()
{
}

coordinate cross(coordinate a, coordinate b)
{
    coordinate c;
    c.x = (a.y*b.z - a.z*b.y);
    c.y = (a.z*b.x - a.x*b.z);
    c.z = (a.x*b.y - a.y*b.x);
    return c;
}

class element
{
public:
    std::vector <coordinate> nodes;
    std::vector <double> plane;
    element();
    element(coordinate a, coordinate b, coordinate c);
    ~element();
};

element::element()
{
    std::vector <coordinate> nodes (3,coordinate());
    std::vector <double> plane (4,0);
    
    this-> nodes = nodes;
    this-> plane = plane;
    

}

element::element(coordinate a, coordinate b, coordinate c)
{
    std::vector <coordinate> nodes {a,b,c};
    std::vector <double> plane (4,0);
    this-> nodes = nodes;
    this-> plane = plane;
}

element::~element()
{
}

class pixel
{
public:
double pixel_val;
    coordinate pixel_centre;
    pixel();
    pixel(coordinate a);
    ~pixel();
};

pixel::pixel()
{
    this-> pixel_val = 0;
    this-> pixel_centre.x = 0;
    this-> pixel_centre.y = 0;
    this-> pixel_centre.z = 0;
}

pixel::pixel(coordinate a)
{
    this-> pixel_val = 0;
    this-> pixel_centre.x = a.x;
    this-> pixel_centre.y = a.y;
    this-> pixel_centre.z = a.z;
}

pixel::~pixel()
{
}

std::vector<std::string> split(const std::string& s)
{
   std::vector<std::string> tokens;
   std::string token;
   tokens.reserve(9);
   std::istringstream tokenStream(s);
   while (getline(tokenStream, token,' '))
   {
      tokens.emplace_back(token);
   }
   return tokens;
}

void DomainParser(std::vector <coordinate> &node_vector,std::vector <element> &element_vector,unsigned int num_of_nodes,unsigned int num_of_elements)
{
    unsigned int counter = 1;
    unsigned int counter1 = 0;
    double a1,a2,a3;
    int ia1,ia2,ia3;
    std::string line;
    std::ifstream infile; 
    infile.open("domain");
    std::vector <std::string> data;
    while (getline(infile,line)) 
    {  
        if (counter == 5)
        {
            data = split(line);
            num_of_nodes = std::stoi(data[0]);
            std::cout << "num_of_nodes = " << num_of_nodes << std::endl;
            node_vector.reserve(num_of_nodes);
        }

        if ( counter > 5 && counter < num_of_nodes+6)
        {
            data = split(line);
            a1 = std::stod(data[1]);
            a2 = std::stod(data[2]);
            a3 = std::stod(data[3]);
            node_vector.emplace_back(coordinate(a1,a2,a3));
        }

        if (counter == num_of_nodes+8)
        {
            data = split(line);
            num_of_elements = std::stoi(data[0]);
            std::cout << "num_of_elements = " << num_of_elements << std::endl;
            element_vector.reserve(num_of_elements);
        }

        if ( counter > num_of_nodes+8 && counter < num_of_nodes+9+num_of_elements)
        {
            data = split(line);
            if (std::stoi(data[1]) == 2)
            {
                ia1 = std::stoi(data[5])-1;
                ia2 = std::stoi(data[6])-1;
                ia3 = std::stoi(data[7])-1;
                element_vector.emplace_back(element(node_vector[ia1],node_vector[ia2],node_vector[ia3]));
                counter1++;
            }    
        }
        counter++;
    }

    element_vector.resize(counter1);
    infile.close();
    //CalcPlane(element_vector);
}

std::vector<std::vector<pixel>> CreateImagePlane(coordinate &cam,double dim1, double dim2, double dim3_coordinate, unsigned int num1, unsigned int num2)
{
    std::vector<std::vector<pixel>> temp_vec (num2,std::vector<pixel>(num1,pixel()));
    double delta1 = dim1/(num1-1);
    double delta2 = dim2/(num2-1);
    for (int i=0;i<num2;i++)
    {
        for (int j=0;j<num2;j++)
        {
            temp_vec[i][j].pixel_centre.x = cam.x + (-dim2*0.5) + ((i)*delta2);
            temp_vec[i][j].pixel_centre.y = cam.y + (-dim1*0.5) + ((j)*delta1);
            temp_vec[i][j].pixel_centre.z = cam.z + dim3_coordinate;
        }
    }

    return temp_vec;
}

void CalcPlane (std::vector <element> &ele_vec)
{
    unsigned int count = ele_vec.size();
    coordinate r1,r2,nor;
    for (int i = 0; i < count; i++)
    {
        //std::cout << "element01 = " << ele_vec[i].nodes[0].x << ele_vec[i].nodes[0].y << ele_vec[i].nodes[0].z << std::endl;
        //std::cout << "element02 = " << ele_vec[i].nodes[1].x << ele_vec[i].nodes[1].y << ele_vec[i].nodes[1].z << std::endl;
        //std::cout << "element03 = " << ele_vec[i].nodes[2].x << ele_vec[i].nodes[2].y << ele_vec[i].nodes[2].z << std::endl; 
        r1 = ele_vec[i].nodes[1] - ele_vec[i].nodes[0];
        r2 = ele_vec[i].nodes[2] - ele_vec[i].nodes[0];
        //std::cout << "r1 = " << r1.x << r1.y << r1.z << std::endl;
        //std::cout << "r2 = " << r2.x << r2.y << r2.z << std::endl;    
        
        nor = cross(r1,r2);
        ele_vec[i].plane[0] = nor.x;
        ele_vec[i].plane[1] = nor.y; 
        ele_vec[i].plane[2] = nor.z;
        ele_vec[i].plane[3] = -(nor.x*ele_vec[i].nodes[0].x + nor.y*ele_vec[i].nodes[0].y + nor.z*ele_vec[i].nodes[0].z);
    }    
}

coordinate CalcInterPoint(element &ele, ray &r)
{
    double t = (ele.plane[0]*r.coeff_x[0] + ele.plane[1]*r.coeff_y[0] + ele.plane[2]*r.coeff_z[0] + ele.plane[3])/
                (ele.plane[0]*r.coeff_x[1] + ele.plane[1]*r.coeff_y[1] + ele.plane[2]*r.coeff_z[1]);
    //std::cout << "t = " << t << std::endl;
    coordinate res;
    res.x = r.coeff_x[0] - t*r.coeff_x[1];
    res.y = r.coeff_y[0] - t*r.coeff_y[1];
    res.z = r.coeff_z[0] - t*r.coeff_z[1];
    return res;
}

bool WithinElement(element &ele, coordinate &p)
{
    using namespace std;
    using namespace Eigen;
    Matrix3f A;
    Vector3f b;
    A << ele.nodes[0].x,ele.nodes[1].x,ele.nodes[2].x,
         ele.nodes[0].y,ele.nodes[1].y,ele.nodes[2].y,
         ele.nodes[0].z,ele.nodes[1].z,ele.nodes[2].z;
    b << p.x,p.y,p.z;
    //cout << "Here is the matrix A:\n" << A << endl;
    //cout << "Here is the vector b:\n" << b << endl;
    Vector3f x = A.colPivHouseholderQr().solve(b);
    //Vector3f x = A.llt().solve(b);
    //cout << "The solution is:\n" << x << endl;
    if (x(0) > 0 && x(1) > 0 && x(2) > 0){
        //cout << "The solution is:\n" << x << endl;
        return true;}
    else{return false;}    
}

void WriteImage(std::vector<std::vector<pixel>> &ip, int k=0, int precesion = 8)
{
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "image-" + s +".dat";
    write.open (path);

    int num1 = ip.size();
    int num2 = ip[0].size();

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << ip[i][j].pixel_val << " ";
                
            }
            else
            {
                write << std::setprecision(precesion) << ip[i][j].pixel_val << "\n";
                
            }
            
        }

    }
    
    write.close();
}

void WriteImagePlane(std::vector<std::vector<pixel>> &ip, int k=0, int precesion = 8)
{
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "data-" + s +".dat";
    write.open (path);

    int num1 = ip.size();
    int num2 = ip[0].size();

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << ip[i][j].pixel_centre.x << " ";
            }
            else
            {
                write << std::setprecision(precesion) << ip[i][j].pixel_centre.x << "\n";
            }
            
        }

    }
    
    write.close();
}

bool check (coordinate lght,coordinate a, coordinate b)
{
    coordinate c = lght - a;
    coordinate d = b - a;
    coordinate e = d/c;
    if (e.x > 0 && e.y > 0 && e.z >0){return true;}
    else{return false;}
}


void render(std::vector <element> &ev,std::vector<std::vector<pixel>> &ip,coordinate &cam,coordinate &light)
{
    
    double light_intensity = 10;
    //std::cout << "cam " << cam.x << cam.y << cam.z << std::endl;
    coordinate inter_secion, surface_point;
    int num_of_int = 0;
    double dis;
    double min_dis;
    int ele_num;
    unsigned int count1 = ip.size();
    unsigned int count2 = ip[0].size();
    std::cout << "count2 = " << count2 << std::endl;
    unsigned int count3 = ev.size();
    unsigned int count4 = ev.size();


    for (int i = 0; i < count1; i++)
    {
        for (int j = 0; j < count2; j++)
        {
            min_dis = 1e10;
            //std::cout << j << std::endl;
            ray r(cam,ip[i][j].pixel_centre);
            
            for (int k = 0; k < count3; k++)
            {
                inter_secion = CalcInterPoint(ev[k],r);
                if(WithinElement(ev[k],inter_secion))
                {
                    //ip[i][j].pixel_val += 1;
                    dis = distance(inter_secion,cam);
                    if(dis<min_dis)
                    {
                        
                        r.tag = 1;
                        min_dis = dis;
                        surface_point = inter_secion;
                        ele_num = k;
                        //std::cout << "dis = " << dis << std::endl;
                        //std::cout << "surface_point = " << surface_point.x << surface_point.y << surface_point.z << std::endl;
                        
                    }
                }
            }
            //std::cout << ele_num << std::endl;
            if(r.tag == 1)   
            {
                ray r1(surface_point,light);
                for (int ii = 0; ii < count4; ii++)
                {
                    if (ii != ele_num)
                    {
                        inter_secion = CalcInterPoint(ev[ii],r1);
                        if(WithinElement(ev[ii],inter_secion))
                        {
                            if(check(light,surface_point,inter_secion))
                            {
                                //std::cout << "sdhf = " << std::endl;
                                //std::cout << "inter_secion = " << inter_secion.x << inter_secion.y << inter_secion.z << std::endl;
                                r1.tag = 2;
                            }
                        }
                    }
                    
                }

                if(r1.tag != 2)
                {
                    //std::cout << j << std::endl;
                    //std::cout << i << j << std::endl;
                    ip[i][j].pixel_val = 10/min_dis/min_dis;
                }
            }   
        }   
    }
    std::cout << "Rendering done.\n" << std::endl;   
}


int main()
{
    std::vector <coordinate> node_vector;
    std::vector <element> element_vector;
    unsigned int num_of_nodes;
    unsigned int num_of_elements;

    std::vector<std::vector<pixel>> image_plane;
    coordinate cam(10,10,0);
    coordinate light(10,10,0);
    double dim1 = 1;
    double dim2 = 1;
    double dim3_coordinate = -5;
    unsigned int num1 = 50;
    unsigned int num2 = 50;


    DomainParser(node_vector,element_vector,num_of_nodes,num_of_elements);
    CalcPlane(element_vector);
    image_plane = CreateImagePlane(cam,dim1,dim2,dim3_coordinate,num1,num2);
    WriteImagePlane(image_plane);
    render(element_vector,image_plane,cam,light);
    WriteImage(image_plane);
    return 0;
}
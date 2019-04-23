#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <random>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <math.h>

class particle;
class species;// public particle;
class mesh;


class particle
{
public:
    std::vector <double> position;
    std::vector <double> velocity;
    std::vector <double> GetPosition(){return position;}
    void SetPosition(std::vector <double> position){position = position;}
    std::vector <double> GetVelocity(){return velocity;}
    void SetVelocity(std::vector <double> velocity){velocity = velocity;}
    particle();
    particle(double, double);
    ~particle();  
};

//Define default Constructer
particle::particle()
{
    std::vector <double> position (2,0.0);
    std::vector <double> velocity (2,0.0);
    this->position = position;
    this->velocity = velocity;
    //std::cout << " A particle is created" << std::endl;
}

//Parameterized Constructor
particle::particle(double charge, double mass)
{
    std::vector <double> position (2,0.0);
    std::vector <double> velocity (2,0.0);
    //std::cout << " A parametrized particle is created" << std::endl;
}

particle::~particle()
{
    //std::cout << " A particle is destroyed";
}

class species: public particle
{
private:
    std::string species_name;
    int species_identifier;
    double charge;
    double mass;
public:
    species();
    species(std::string, int, double, double);
    ~species();
    double GetCharge(){return charge;}
    void SetCharge(double charge){this->charge=charge;}
    double GetMass(){return mass;}
    void SetMass(double mass){this->mass = mass;}
    int GetID(){return species_identifier;}
    void SetID(int id){this->species_identifier=id;}
    void print_properties();
};

species::species() : particle()
{
    this->species_name = "no name";
    this->species_identifier = 0;
    this->charge = 0;
    this->mass = 1;
    //std::cout << " A species particle is created" << std::endl;
}


species::species(std::string species_name , int species_identifier, double charge,
double mass) : particle()
{
    this->species_name = species_name;
    this->species_identifier = species_identifier;
    this->charge = charge;
    this->mass = mass;
    //std::cout << " A parametrized species is created" << std::endl;
}

species::~species()
{
    //std::cout << " A species particle is destroyed" << std::endl;
}

void species::print_properties()
{
    using std::cout;
    using std::endl;
    cout << "charge = " << this->charge << endl;
    cout << "mass = " << this->mass << endl;
    cout << "position_x = " << position[0] <<endl;
    cout << "velocity_x = " << velocity[0] <<endl;
}

class mesh
{
//private:
//    int num_dim1;
//    int num_dim2;
//    double h_dim1;
//    double h_dim2;
//    
public:
    double charge;
    double nodal_volume;
    double potential;
    double density_0;
    double density_1;
    std::vector <double> electric_feild;
    std::vector <double> magnetic_feild;
    mesh();
    ~mesh();
};

mesh::mesh()
{
    std::vector <double> electric_feild (2,0);
    std::vector <double> magnetic_feild (2,0);
    this->charge = 0;
    this->nodal_volume = 0;
    this->potential = 0;
    this->density_0 = 0;
    this->density_1 = 0;
    this->electric_feild = electric_feild;
    this->magnetic_feild = magnetic_feild;
}

mesh::~mesh()
{
}

struct species_data
{
    std::string species_name;
    int species_id;
    double species_charge;
    double _species_mass;
};

std::vector <species> create_species (struct species_data data, double num_of_particles = 0)
{
    //std::cout << "done" << std::endl;
    std::vector <species> temp_array (num_of_particles,species(data.species_name,data.species_id,data.species_charge,data._species_mass));
    return temp_array;
}

std::vector <species> load_species (struct species_data data,double min_dim1, double max_dim1, double min_dim2, double max_dim2, double density)
{
    double delta_dim1 = (max_dim1-min_dim1);
    double delta_dim2 = (max_dim2-min_dim2);
    int num_of_particles = delta_dim1*((max_dim2*max_dim2)-(min_dim2*min_dim2))*density;
    //std::cout << num_of_particles << std::endl;
    int N = 50;
    double dy = (max_dim2-min_dim2)/N;
    std::vector <species> temp_array1 = create_species(data,num_of_particles);
    int counter1 = 0;
    int counter2 = 0;
    for (int i = 0; i<N; i++)
    {
        double r1 = (i*dy) + min_dim2;
        double r2 = ((i+1)*dy) + min_dim2;
        int num = ((r2*r2) - (r1*r1))*delta_dim1*density;
        counter2 += num;
        for (int j = counter1; j<counter2; j++)
        {
            temp_array1[j].position[0] =  min_dim1 + ((rand() / (RAND_MAX + 1.))*delta_dim1);
            temp_array1[j].position[1] = r1 + ((rand() / (RAND_MAX + 1.))*dy);
            if (temp_array1[j].position[1] == -1)
            {
                //std::cout << i << std::endl;
            }
        }
        counter1 = counter2;   
    }
    //std::cout << counter2 << std::endl;
    if (counter2 < num_of_particles)
    {
        int num = num_of_particles - counter2;
        counter2 += num;
        double r1 = (0.65*delta_dim2) + min_dim2;
        double r2 = (0.99*delta_dim2) + min_dim2;
        for (int i = counter1; i < counter2; i++)
        {
            temp_array1[i].position[0] =  min_dim1 + ((rand() / (RAND_MAX + 1.))*delta_dim1);
            temp_array1[i].position[1] = r1 + ((rand() / (RAND_MAX + 1.))*(r2-r1)); 
        }
    }
    //std::cout << counter2 << std::endl;
    
    return temp_array1;
}

void set_thermal_velocity(std::vector <species> &species_array, double ThermalVelocity = 0, int seed = 4444)
{   
    std::default_random_engine e(seed);
    std::normal_distribution <double> distrN (1.0,1.0);
    int N = species_array.size();
    for (int i=0; i<N; i++)
    {
        species_array[i].velocity[0] = distrN(e)*ThermalVelocity;
        species_array[i].velocity[1] = distrN(e)*ThermalVelocity;
    }
}

void set_beam_velocity(std::vector <species> &species_array, std::vector <double> &BeamVelocity)
{   
    int N = species_array.size();
    for (int i=0; i<N; i++)
    {
        species_array[i].velocity[0] += BeamVelocity[0];
        species_array[i].velocity[1] += BeamVelocity[1];
    }
}

std::vector< std::vector<mesh> > CreateMesh(double min_dim1, double max_dim1, double min_dim2, double max_dim2, int num1, int num2)
{
    std::vector<std::vector<mesh>> temp_array(num1+2, std::vector<mesh>(num2+2,mesh()));
    return temp_array;
}

void calc_nodal_volume(std::vector< std::vector<mesh> > &mesh_array, double min_dim2, double h_dim1, double h_dim2)
{
    int num1 = mesh_array.size() - 2;
    int num2 = mesh_array[0].size() - 2;

    for (int i=1;i<num1+1;i++)
    {
        for (int j=1;j<num2+1;j++)
        {
            double j_min = (j-1) - 0.5;
            double j_max = (j-1) + 0.5;
            if (j_min<0){j_min=0;}
            if (j_max>num2){j_max=num2;}
            double factor;
            if (i == 1 || i == num1){factor = 0.5;}
            else{factor = 1;}
            double r_min = min_dim2 + (h_dim2*j_min);
            double r_max = min_dim2 + (h_dim2*j_max);
            mesh_array[i][j].nodal_volume = factor*h_dim1*((r_max*r_max)-(r_min*r_min));
        }
    }
    
}

void reset_charge(std::vector< std::vector<mesh>> &mesh_arr)
{
    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            mesh_arr[i][j].charge = 0;
            mesh_arr[i][j].potential = 0;
        }
        
    }
}


void reset_density(std::vector< std::vector<mesh>> &mesh_arr)
{
    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            mesh_arr[i][j].density_0 = 0;
            mesh_arr[i][j].density_1 = 0;
        }
        
    }
}

void charge_weighting(std::vector< std::vector<mesh>> &mesh_arr, std::vector<species> &species_arr,
double min_dim1, double min_dim2, double h_dim1, double h_dim2)
{
    
    int n = species_arr.size();
    double particle_charge = species_arr[0].GetCharge();
    double px,py,hx,hy;
    double a,b,c,d;
    int ix0,iy0,ix1,iy1;
    for (int i=0; i<n; i++)
    {
        //std::cout << "i = "<< i << std::endl;
        px = species_arr[i].position[0] - min_dim1;
        py = species_arr[i].position[1] - min_dim2;
        ix0 = px/h_dim1;
        iy0 = py/h_dim2;
        ix1 = ix0 + 1;
        iy1 = iy0 + 1;
        hx = ((px/h_dim1)-ix0);
        hy = ((py/h_dim2)-iy0);
        a = (1-hx)*(1-hy);
        b = hx*(1-hy);
        c = (1-hx)*hy;
        d = hx*hy;

        mesh_arr[ix0+1][iy0+1].charge += a*particle_charge;
        mesh_arr[ix1+1][iy0+1].charge += b*particle_charge;
        mesh_arr[ix0+1][iy1+1].charge += c*particle_charge;
        mesh_arr[ix1+1][iy1+1].charge += d*particle_charge;
        
    }
    
}

void density0_weighting(std::vector< std::vector<mesh>> &mesh_arr, std::vector<species> &species_arr,
double min_dim1, double min_dim2, double h_dim1, double h_dim2)
{
    
    int n = species_arr.size();
    //double particle_charge = species_arr[0].GetCharge();
    for (int i=0; i<n; i++)
    {
        //std::cout << "i = "<< i << std::endl;
        double px = species_arr[i].position[0] - min_dim1;
        double py = species_arr[i].position[1] - min_dim2;
        int ix0 = px/h_dim1;
        int iy0 = py/h_dim2;
        int ix1 = ix0 + 1;
        int iy1 = iy0 + 1;
        double hx = ((px/h_dim1)-ix0);
        double hy = ((py/h_dim2)-iy0);
        mesh_arr[ix0+1][iy0+1].density_0 += (1-hx)*(1-hy);
        mesh_arr[ix1+1][iy0+1].density_0 += hx*(1-hy);
        mesh_arr[ix0+1][iy1+1].density_0 += (1-hx)*hy;
        mesh_arr[ix1+1][iy1+1].density_0 += hx*hy;
        
    }
    
}

void density1_weighting(std::vector< std::vector<mesh>> &mesh_arr, std::vector<species> &species_arr,
double min_dim1, double min_dim2, double h_dim1, double h_dim2)
{
    
    int n = species_arr.size();
    //double particle_charge = species_arr[0].GetCharge();
    for (int i=0; i<n; i++)
    {
        //std::cout << "i = "<< i << std::endl;
        double px = species_arr[i].position[0] - min_dim1;
        double py = species_arr[i].position[1] - min_dim2;
        int ix0 = px/h_dim1;
        int iy0 = py/h_dim2;
        int ix1 = ix0 + 1;
        int iy1 = iy0 + 1;
        double hx = ((px/h_dim1)-ix0);
        double hy = ((py/h_dim2)-iy0);
        mesh_arr[ix0+1][iy0+1].density_1 += (1-hx)*(1-hy);
        mesh_arr[ix1+1][iy0+1].density_1 += hx*(1-hy);
        mesh_arr[ix0+1][iy1+1].density_1 += (1-hx)*hy;
        mesh_arr[ix1+1][iy1+1].density_1 += hx*hy;
        
    }
    
}


void possion_solve(std::vector< std::vector<mesh>> &mesh_arr,double h_dim1, double h_dim2, double epsilon, double w)
{
    int num1 = mesh_arr.size() - 2;
    int num2 = mesh_arr[0].size() - 2;

    for (int j=1;j<num2+1;j++)
    {
        mesh_arr[1][j].potential =  mesh_arr[2][j].potential;//200;       // SET LEFT BOUNDARY CONDITIONS
    }
    
    for (int j=1;j<num2+1;j++)
    {
        mesh_arr[num1][j].potential = mesh_arr[num1-1][j].potential;//-200;    // SET RIGHT BOUNDARY CONDITIONS
    }

    for (int i=1;i<num1+1;i++)
    {
        mesh_arr[i][1].potential = mesh_arr[i][2].potential;
    }

    for (int i=1;i<num1+1;i++)
    {
        mesh_arr[i][num2].potential = 0;//mesh_arr[i][num2-1].potential;
    }
    
    for (int i=2;i<num1;i++)
    {    
        for (int j=2;j<num2;j++)
        {
            if ((i+j)%2 == 0)
            {
                double prev = mesh_arr[i][j].potential;
                double t1 = mesh_arr[i][j].charge/mesh_arr[i][j].nodal_volume/epsilon;
                double t2 = (mesh_arr[i][j+1].potential + mesh_arr[i][j-1].potential)/(h_dim2*h_dim2);
                double t3 = (mesh_arr[i][j+1].potential - mesh_arr[i][j-1].potential)/(2.0*h_dim2*h_dim2*(j-1));
                double t4 = (mesh_arr[i+1][j].potential + mesh_arr[i-1][j].potential)/(h_dim1*h_dim1);
                double t5 = (2.0/(h_dim2*h_dim2)) + (2.0/(h_dim1*h_dim1));
                mesh_arr[i][j].potential = prev + w*(((t1+t2+t4+t3)/t5)-prev);
                
            }
        }
        
    }
    
    for (int i=2;i<num1;i++)
    {    
        for (int j=2;j<num2;j++)
        {
            if ((i+j)%2 != 0)
            {
                double prev = mesh_arr[i][j].potential;
                double t1 = mesh_arr[i][j].charge/mesh_arr[i][j].nodal_volume/epsilon;
                double t2 = (mesh_arr[i][j+1].potential + mesh_arr[i][j-1].potential)/(h_dim2*h_dim2);
                double t3 = (mesh_arr[i][j+1].potential - mesh_arr[i][j-1].potential)/(2.0*h_dim2*h_dim2*(j-1));
                double t4 = (mesh_arr[i+1][j].potential + mesh_arr[i-1][j].potential)/(h_dim1*h_dim1);
                double t5 = (2.0/(h_dim2*h_dim2)) + (2.0/(h_dim1*h_dim1));
                mesh_arr[i][j].potential = prev + w*(((t1+t2+t4+t3)/t5)-prev);
            }
        }
        
    }

}

void update_ghost_cells(std::vector< std::vector<mesh>> &mesh_arr)
{
    int num1 = mesh_arr.size() - 2;
    int num2 = mesh_arr[0].size() - 2;

    for (int i=1;i<num1+1;i++)
    {
        mesh_arr[i][0].potential = 2.0*mesh_arr[i][1].potential - mesh_arr[i][2].potential;
    }

    for (int i=1;i<num1+1;i++)
    {
        mesh_arr[i][num2+1].potential = (2.0*mesh_arr[i][num2].potential) - (mesh_arr[i][num2-1].potential);
    }

    for (int j=0;j<num2+2;j++)
    {
        mesh_arr[num1+1][j].potential = 2.0*mesh_arr[num1][j].potential - mesh_arr[num1-1][j].potential;
    }

    for (int j=0;j<num2+2;j++)
    {
        mesh_arr[0][j].potential = 2.0*mesh_arr[1][j].potential - mesh_arr[2][j].potential;
    }
}

void calculate_electric_feild(std::vector< std::vector<mesh>> &mesh_arr, double h_dim1, double h_dim2)
{
    int num1 = mesh_arr.size() - 2;
    int num2 = mesh_arr[0].size() - 2;

    for (int i=1;i<num1+1;i++)
    {    
        for (int j=1;j<num2+1;j++)
        {
            mesh_arr[i][j].electric_feild[0] = (mesh_arr[i-1][j].potential - mesh_arr[i+1][j].potential)/(2*h_dim1);
            mesh_arr[i][j].electric_feild[1] = (mesh_arr[i][j-1].potential - mesh_arr[i][j+1].potential)/(2*h_dim2);
        }
    }
}

void interpolate_feild(double &px, double &py, double &ex, double &ey,std::vector< std::vector<mesh>> &mesh_arr,
double min_dim1, double min_dim2,double h_dim1, double h_dim2)
{
    int ix0,iy0,ix1,iy1;
    double hx,hy,a,b,c,d,val0,val1,val2,val3;

    ix0 = px/h_dim1;
    iy0 = py/h_dim2;
    ix1 = ix0 + 1;
    iy1 = iy0 + 1;
    hx = ((px/h_dim1)-ix0);
    hy = ((py/h_dim2)-iy0);
    a = (1-hx)*(1-hy);
    b = hx*(1-hy);
    c = (1-hx)*hy;
    d = hx*hy;
    
    val0 = a*mesh_arr[ix0][iy0].electric_feild[0];
    val1 = b*mesh_arr[ix1][iy0].electric_feild[0];
    val2 = c*mesh_arr[ix0][iy1].electric_feild[0];
    val3 = d*mesh_arr[ix1][iy1].electric_feild[0];
    ex = (val0+val1+val2+val3);
    val0 = a*mesh_arr[ix0][iy0].electric_feild[1];
    val1 = b*mesh_arr[ix1][iy0].electric_feild[1];
    val2 = c*mesh_arr[ix0][iy1].electric_feild[1];
    val3 = d*mesh_arr[ix1][iy1].electric_feild[1];
    ey = (val0+val1+val2+val3);
}

void boris_pusher(std::vector <species> &species_arr,std::vector< std::vector<mesh>> &mesh_arr, double min_dim1, double min_dim2,
double h_dim1, double h_dim2 , double time_step)
{
    int num = species_arr.size();
    double particle_chrg,particle_mass;
    double px,py,ux,vy,E0,E1;
    //double x1,y1,u1,v1,x2,y2,u2,v2,x3,y3,u3,v3,x4,y4,u4,v4;
    //double dx1,dy1,du1,dv1,dx2,dy2,du2,dv2,dx3,dy3,du3,dv3,dx4,dy4,du4,dv4;
    int ix0,iy0,ix1,iy1;
    double hx,hy,a,b,c,d,val0,val1,val2,val3;

    for (int i=0; i<num; i++)
    {   
        
        particle_chrg = species_arr[i].GetCharge();
        particle_mass = species_arr[i].GetMass();
        px = species_arr[i].position[0];
        py = species_arr[i].position[1];
        ux = species_arr[i].velocity[0];
        vy = species_arr[i].velocity[1];
        //interpolate_feild(px,py,E0,E1,mesh_arr,min_dim1,min_dim2,h_dim1,h_dim2);
        ix0 = px/h_dim1;
        iy0 = py/h_dim2;
        ix1 = ix0 + 1;
        iy1 = iy0 + 1;
        hx = ((px/h_dim1)-ix0);
        hy = ((py/h_dim2)-iy0);
        a = (1-hx)*(1-hy);
        b = hx*(1-hy);
        c = (1-hx)*hy;
        d = hx*hy;
    
        val0 = a*mesh_arr[ix0][iy0].electric_feild[0];
        val1 = b*mesh_arr[ix1][iy0].electric_feild[0];
        val2 = c*mesh_arr[ix0][iy1].electric_feild[0];
        val3 = d*mesh_arr[ix1][iy1].electric_feild[0];
        E0 = (val0+val1+val2+val3);
        val0 = a*mesh_arr[ix0][iy0].electric_feild[1];
        val1 = b*mesh_arr[ix1][iy0].electric_feild[1];
        val2 = c*mesh_arr[ix0][iy1].electric_feild[1];
        val3 = d*mesh_arr[ix1][iy1].electric_feild[1];
        E1 = (val0+val1+val2+val3);
        
        species_arr[i].velocity[0] += particle_chrg*E0*time_step/particle_mass;
        species_arr[i].velocity[1] += particle_chrg*E1*time_step/particle_mass;
        species_arr[i].position[0] += species_arr[i].velocity[0]*time_step;
        species_arr[i].position[1] += species_arr[i].velocity[1]*time_step;
        
        /*
        dx1 = ux*time_step;
        dy1 = vy*time_step;
        du1 = particle_chrg*E0*time_step/particle_mass;
        dv1 = particle_chrg*E1*time_step/particle_mass;

        x1 = px + dx1*0.5;
        y1 = py + dy1*0.5;
        u1 = ux + du1*0.5;
        v1 = vy + dv1*0.5;
        interpolate_feild(x1,y1,E0,E1,mesh_arr,min_dim1,min_dim2,h_dim1,h_dim2);
        dx2 = u1*time_step;
        dy2 = v1*time_step;
        du2 = particle_chrg*E0*time_step/particle_mass;
        dv2 = particle_chrg*E1*time_step/particle_mass;
        
        x2 = px + dx2*0.5;
        y2 = py + dy2*0.5;
        u2 = ux + du2*0.5;
        v2 = vy + dv2*0.5;
        interpolate_feild(x2,y2,E0,E1,mesh_arr,min_dim1,min_dim2,h_dim1,h_dim2);

        dx3 = u2*time_step;
        dy3 = v2*time_step;
        du3 = particle_chrg*E0*time_step/particle_mass;
        dv3 = particle_chrg*E1*time_step/particle_mass;
        x3 = px + dx3;
        y3 = py + dy3;
        u3 = ux + du3;
        v3 = vy + dv3;
        interpolate_feild(x3,y3,E0,E1,mesh_arr,min_dim1,min_dim2,h_dim1,h_dim2);

        dx4 = u3*time_step;
        dy4 = v3*time_step;
        du4 = particle_chrg*E0*time_step/particle_mass;
        dv4 = particle_chrg*E1*time_step/particle_mass;

        species_arr[i].position[0] += ((dx1+2*dx2+2*dx3+dx4)/6.0);
        species_arr[i].position[1] += ((dy1+2*dy2+2*dy3+dy4)/6.0);
        species_arr[i].velocity[0] += ((du1+(2*du2)+(2*du3)+du4)/6.0);
        species_arr[i].velocity[1] += ((dv1+2*dv2+2*dv3+dv4)/6.0);
        //std::cout << du1 << du2 << du3 << du4 << species_arr[i].velocity[0] << std::endl;
        

        x2 = px + 0.5*ux*time_step;
        y2 = py + 0.5*vy*time_step;
        u2 = ux + 0.5*E0*time_step*particle_chrg/particle_mass;
        v2 = vy + 0.5*E1*time_step*particle_chrg/particle_mass;
        interpolate_feild(x2,y2,E0,E1,mesh_arr,min_dim1,min_dim2,h_dim1,h_dim2);
        Ex1 = E0;
        Ey1 = E1;


        x3 = px + 0.5*u2*time_step;
        y3 = py + 0.5*v2*time_step;
        u3 = ux + 0.5*E0*time_step*particle_chrg/particle_mass;
        v3 = vy + 0.5*E1*time_step*particle_chrg/particle_mass;
        interpolate_feild(x3,y3,E0,E1,mesh_arr,min_dim1,min_dim2,h_dim1,h_dim2);
        Ex2 = E0;
        Ey2 = E1;

        x4 = px + u3*time_step;
        y4 = py + v3*time_step;
        u4 = ux + E0*time_step*particle_chrg/particle_mass;
        v4 = vy + E1*time_step*particle_chrg/particle_mass;
        interpolate_feild(x3,y3,E0,E1,mesh_arr,min_dim1,min_dim2,h_dim1,h_dim2);
        Ex3 = E0;
        Ey3 = E1;

        species_arr[i].position[0] += ((ux+2*u2+2*u3+u4)/6.0)*time_step;
        species_arr[i].position[1] += ((vy+2*v2+2*v3+v4)/6.0)*time_step;
        species_arr[i].velocity[0] += ((Ex0+2*Ex1+2*Ex2+Ex3)/6.0)*time_step;
        species_arr[i].velocity[1] += ((Ey0+2*Ey1+2*Ey2+Ey3)/6.0)*time_step;*/
    }
}

void merge(std::vector <species> &species_arr1,std::vector <species> &species_arr2)
{
    int num1 = species_arr1.size();
    int num2 = species_arr2.size();
    int num3 = num1 + num2;
    //std::cout << num1 << std::endl;
    species_arr1.reserve(num3);
    species_arr1.insert(species_arr1.end(), species_arr2.begin(), species_arr2.end());
    //std::cout << "c = " << species_arr1.size() << std::endl;
}

void source(std::vector <species> &species_arr,struct species_data data,double min_dim1, double max_dim1, double min_dim2, double max_dim2, double density, double thermalvelocity, std::vector <double> &BeamVelocity)
{
    std::vector <species> arr = load_species(data,min_dim1, max_dim1, min_dim2, max_dim2, density);
    //std::cout << "a = " << arr.size() << std::endl;
    set_thermal_velocity(arr,thermalvelocity);
    set_beam_velocity(arr,BeamVelocity);
    merge(species_arr,arr);
    //std::cout << "b = " << species_arr.size() << std::endl;
}

void sink_r(std::vector <species> &species_arr, double coordinate, const int p)
{
    for (int i=0; i<species_arr.size();i++)
    {
        if (species_arr[i].position[p] >= coordinate)
        {
            species_arr.erase(species_arr.begin()+i);
            --i;

        }
    }
}

void sink_l(std::vector <species> &species_arr, double coordinate, int p)
{
    int num = species_arr.size();
    for (int i=0; i<species_arr.size();i++)
    {
        if (species_arr[i].position[p] <= coordinate)
        {
            species_arr.erase(species_arr.begin()+i);
            --i;
        }
    }
}

void mirror_u(std::vector <species> &species_arr, double coordinate, int p)
{
    int num = species_arr.size();
    for (int i=0; i<num;i++)
    {
        if (species_arr[i].position[p] >= coordinate)
        {
            //std::cout << "before" << species_arr[i].position[p] << "\n";
            species_arr[i].position[p] = coordinate - abs(species_arr[i].position[p]-coordinate);
            //std::cout << "during" << coordinate - abs(species_arr[i].position[p]-coordinate) << "\n";
            species_arr[i].velocity[p] = -species_arr[i].velocity[p];
            //std::cout << "after" << species_arr[i].position[p] << "\n";
        }
    }
}

void mirror_b(std::vector <species> &species_arr, double coordinate, int p)
{
    int num = species_arr.size();
    for (int i=0; i<num;i++)
    {
        
        if (species_arr[i].position[p] <= coordinate)
        {
            species_arr[i].position[p] = coordinate + abs(species_arr[i].position[p]-coordinate);
            species_arr[i].velocity[p] = -species_arr[i].velocity[p];
        }
    }
}

void calc_charge_on_wall_l(std::vector <species> &species_arr,double &chrg,double coordinate, int p)
{
    int num = species_arr.size();
    for (int i=0; i<num;i++)
    {
        if (species_arr[i].position[p] <= coordinate)
        {
            chrg += species_arr[i].GetCharge();
        }
    }

    sink_l(species_arr,coordinate,p);
}

void calc_charge_on_wall_u(std::vector <species> &species_arr,double &chrg,double coordinate, int p)
{
    int num = species_arr.size();
    for (int i=0; i<num;i++)
    {
        if (species_arr[i].position[p] >= coordinate)
        {
            chrg += species_arr[i].GetCharge();
        }
    }

    sink_r(species_arr,coordinate,p);
}

void charge_to_mesh_l(std::vector< std::vector<mesh>> &mesh_arr,double &chrg, double coordinate, double min_dim2 , double h_dim2)
{
    int j = (coordinate - min_dim2)/h_dim2;
    int num1 = mesh_arr.size();

    for (int i=1;i<num1-1;i++)
    {
        mesh_arr[i][j+1].charge += chrg;
    }
}

void write_output_particles(std::vector <species> &species_arr1,std::vector <species> &species_arr2, int k, int precesion = 8)
{
    /*system("mkdir Output");
    system("mkdir Output/Particles");
    system("mkdir Output/Charge");
    system("mkdir Output/ElectricFeild");
    system("mkdir Output/ElectricFeild/dim1");
    system("mkdir Output/ElectricFeild/dim2");
    system("mkdir Output/MagneticFeild");
    system("mkdir Output/MagneticFeild/dim1");
    system("mkdir Output/MagneticFeild/dim2");
    system("mkdir Output/PotentialFeild");
    */
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "Output/Particles/phase-" + s +".dat";
    write.open (path);
    int num1 = species_arr1.size();
    int num2 = species_arr2.size();
    //std::cout << species_arr1[1].velocity[0] << std::endl;
    for (int i=0; i<num1+num2;i++)
    {
        if (i<num1)
        {
            write << std::setprecision(precesion) << species_arr1[i].GetID() << " " << species_arr1[i].GetCharge() << " " << species_arr1[i].GetMass() << " " << species_arr1[i].position[0] << " " << species_arr1[i].position[1] << " " << species_arr1[i].velocity[0] << " " << species_arr1[i].velocity[1] << "\n";
        }
        else
        {
            write << std::setprecision(precesion) << species_arr2[i-num1].GetID() << " " << species_arr2[i-num1].GetCharge() << " " << species_arr2[i-num1].GetMass() << " " << species_arr2[i-num1].position[0] << " " << species_arr2[i-num1].position[1] << " " << species_arr2[i-num1].velocity[0] << " "<< species_arr2[i-num1].velocity[1] << "\n";
        }
        
    }
    write.close();
}

void write_output_charge(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "Output/Charge/charge-" + s +".dat";
    write.open (path);

    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].charge << " ";
            }
            else
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].charge << "\n";
            }
            
        }

    }
    
    write.close();
}

void write_output_density0(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "Output/Density/0/tag0-" + s +".dat";
    write.open (path);

    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();
    double value;
    for (int i=1;i<num1-1;i++)
    {    
        for (int j=1;j<num2-1;j++)
        {
            if (j != num2-2)
            {
                value = mesh_arr[i][j].density_0/mesh_arr[i][j].nodal_volume;
                write << std::setprecision(precesion) << value << " ";
            }
            else
            {
                value = mesh_arr[i][j].density_0/mesh_arr[i][j].nodal_volume;
                write << std::setprecision(precesion) << value << "\n";
            }
            
        }

    }
    
    write.close();
}

void write_output_density1(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "Output/Density/1/tag1-" + s +".dat";
    write.open (path);

    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();
    double value;

    for (int i=1;i<num1-1;i++)
    {    
        for (int j=1;j<num2-1;j++)
        {
            if (j != num2-2)
            {
                value = mesh_arr[i][j].density_1/mesh_arr[i][j].nodal_volume;
                write << std::setprecision(precesion) << value << " ";
            }
            else
            {
                value = mesh_arr[i][j].density_1/mesh_arr[i][j].nodal_volume;
                write << std::setprecision(precesion) << value << "\n";
            }
            
        }

    }
    
    write.close();
}

void write_output_nodal_volume(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "Output/NodalVolume/nv-" + s +".dat";
    write.open (path);

    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].nodal_volume << " ";
            }
            else
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].nodal_volume << "\n";
            }
            
        }

    }
    
    write.close();
}

void write_output_potential(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    std::ofstream write;
    std::string s = std::to_string(k);
    std::string path = "Output/PotentialFeild/potential-" + s +".dat";
    write.open (path);

    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].potential << " ";
            }
            else
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].potential << "\n";
            }
            
        }

    }
    
    write.close();
}

void write_output_electricfeild(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    std::ofstream write;
    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();

    std::string s = std::to_string(k);
    std::string path = "Output/ElectricFeild/dim1/elec-dim1-" + s +".dat";
    write.open (path);

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].electric_feild[0] << " ";
            }
            else
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].electric_feild[0] << "\n";
            }
            
        }

    }
    
    write.close();

    path = "Output/ElectricFeild/dim2/elec-dim2-" + s +".dat";
    write.open (path);

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].electric_feild[1] << " ";
            }
            else
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].electric_feild[1] << "\n";
            }
            
        }

    }
    
    write.close();
}

void write_output_magneticfeild(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    std::ofstream write;
    int num1 = mesh_arr.size();
    int num2 = mesh_arr[0].size();

    std::string s = std::to_string(k);
    std::string path = "Output/MagneticFeild/dim1/mag-dim1-" + s +".dat";
    write.open (path);

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].magnetic_feild[0] << " ";
            }
            else
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].magnetic_feild[0] << "\n";
            }
            
        }

    }
    
    write.close();

    path = "Output/MagneticFeild/dim2/mag-dim2-" + s +".dat";
    write.open (path);

    for (int i=0;i<num1;i++)
    {    
        for (int j=0;j<num2;j++)
        {
            if (j != num2-1)
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].magnetic_feild[1] << " ";
            }
            else
            {
                write << std::setprecision(precesion) << mesh_arr[i][j].magnetic_feild[1] << "\n";
            }
            
        }

    }
    
    write.close();
}

void write_output_all_mesh(std::vector< std::vector<mesh>> &mesh_arr, int k, int precesion = 8)
{
    write_output_charge(mesh_arr,k,precesion);
    write_output_potential(mesh_arr,k,precesion);
    write_output_electricfeild(mesh_arr,k,precesion);
    //write_output_magneticfeild(mesh_arr,k,precesion);
    write_output_density0(mesh_arr,k,precesion);
    write_output_density1(mesh_arr,k,precesion);
}

int main()
{
    time_t start, end;
    double charge_top = 0, charge_bottom = 0;
    double elec_mass = 1;
    double elec_chrg = -1;
    double ion_mass = 100;
    double ion_chrg = 1;
    int elec_tag = 0;
    int ion_tag = 1;
    
    species_data species1 = {"electron",0,-1,1};
    species_data species2 = {"ion",1,1,100};

    double density = 1000000;
    double epsilon = 1000;
    const int k_b = 1.38064852;
    double temperture = 100;
    double thermal_vel_elec = 1.7*sqrt(k_b*temperture/elec_mass);
    double thermal_vel_ion = 1.7*sqrt(k_b*temperture/ion_mass);
    double plasma_frequency_elec = sqrt((density*(elec_chrg*elec_chrg))/(elec_mass*epsilon));
    double plasma_frequency_ion = sqrt((density*(ion_chrg*ion_chrg))/(ion_mass*epsilon));
    double debye_length_elec = thermal_vel_elec/plasma_frequency_elec;
    double debye_length_ion = thermal_vel_ion/plasma_frequency_ion;
    std::cout << " plasma_frequency_elec = " << plasma_frequency_elec << ", debye_length_elec = " << debye_length_elec << std::endl;
    std::cout << " plasma_frequency_ion = " << plasma_frequency_ion << ", debye_length_ion = " << debye_length_ion << std::endl;

    double dt = 1.0/(plasma_frequency_elec*100);
    
    double dz = debye_length_elec/50;
    double dr = debye_length_elec/50;

    double z0 = 0;
    double z1 = 2.0*debye_length_ion;
    double r0 = 0;
    double r1 = 1.0*debye_length_ion;
    double wall_rb = dr;
    double wall_ru = r1-dr;

    double tf = 10*(z1-z0)/thermal_vel_elec;

    int num1 = (z1-z0)/dz + 1;
    int num2 = (r1-r0)/dr + 1;
    std::cout << num1 << " " << num2 << std::endl;
    std::cout << "timestep = " << dt << std::endl;

    double w = 1.4;
    
    std::vector<std::vector<mesh>> grid = CreateMesh(z0,z1,r1,r1,num1,num2);
    
    std::vector <species> electrons = load_species(species1,z0,z1,wall_rb,r1,density);
    std::cout << "done" << std::endl;
    std::vector <species> ions = load_species(species2,z0,z1,wall_rb,r1,density);
    set_thermal_velocity(electrons,thermal_vel_elec);
    set_thermal_velocity(ions,thermal_vel_ion);
    std::cout << "num 0f elec = " << electrons.size() << std::endl;

    calc_nodal_volume(grid,r0,dz,dr);
    write_output_nodal_volume(grid,0);

    start = clock();
    int N = 1;//tf/dt;
    std::cout << "no. of timesteps = "<< N << std::endl;
    for (int k=0;k<N;k++)
    {
        calc_charge_on_wall_l(electrons,charge_bottom,wall_rb,1);
        //std::cout << "cherg = "<< charge_bottom << std::endl;
        reset_charge(grid);
        charge_weighting(grid,electrons,z0,r0,dz,dr);
        charge_weighting(grid,ions,z0,r0,dz,dr);
        charge_to_mesh_l(grid,charge_bottom,wall_rb,r0,dr);
        for (int i=0;i<10000;i++)
        {
            possion_solve(grid,dz,dr,epsilon,w);
        }
        update_ghost_cells(grid);
        calculate_electric_feild(grid,dz,dr);
        boris_pusher(electrons,grid,z0,r0,dz,dr,dt);
        boris_pusher(ions,grid,z0,r0,dz,dr,dt);
        mirror_u(electrons,r1,1);
        mirror_b(ions,wall_rb,1);
        mirror_u(ions,r1,1);
        mirror_b(electrons,z0,0);
        mirror_u(electrons,z1,0);
        mirror_b(ions,z0,0);
        mirror_u(ions,z1,0);
        if (k%20 == 0){
            reset_density(grid);
            write_output_particles(electrons,ions,k+1);
            density0_weighting(grid,electrons,z0,r0,dz,dr);
            density1_weighting(grid,ions,z0,r0,dz,dr);
            write_output_all_mesh(grid,k+1);
            std::cout << "k = " << k << std::endl;
        }
        if (k == N-1){
            write_output_particles(electrons,ions,k+1);
            reset_density(grid);
            density0_weighting(grid,electrons,z0,r0,dz,dr);
            density1_weighting(grid,ions,z0,r0,dz,dr);
            write_output_all_mesh(grid,k+1);
            std::cout << "k = " << k << std::endl;
        }
                    
    }
    end = clock();
    double time_taken = double(end - start);
    std::cout << "runtime = " << time_taken / CLOCKS_PER_SEC << "sec" << "\n";
    return 0;
}
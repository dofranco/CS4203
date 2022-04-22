/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: Daniel Franco dofranco@usc.edu
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>

#include <math.h>
#include <glm/glm.hpp>
#include <cmath>
#include <iostream>
#include <vector>

//using namespace glm;

#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

#define M_PI acos(-1.0)

char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

glm::dvec3 zero_vector = { 0.0, 0.0, 0.0 };

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};

struct Triangle
{
  Vertex v[3];

  glm::dvec3 tri_normal;
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;

  glm::dvec3 GetNormal(const glm::dvec3& intersect_point)
  {
      return (glm::dvec3(intersect_point.x - position[0], intersect_point.y - position[1], intersect_point.z - position[2]) / radius);
  }
};

struct Light
{
  double position[3];
  double color[3];
};

struct Ray {

    //default constructor
    Ray() : ray_color(zero_vector),
        tri_intersect(nullptr),
        sph_intersect(nullptr),
        intersect_point( -1, -1, -1 ),
        t_val(0) 
    {};
    
    //constructor
    Ray(glm::dvec3 in_orig, glm::dvec3 in_dir)
    {
        ray_origin = in_orig;
        ray_direction = glm::normalize(in_dir);
        ray_color = zero_vector;
        tri_intersect = nullptr;
        sph_intersect = nullptr;
        intersect_point = { -1, -1, -1 };
    }

    //setter
    void Set_Ray(glm::dvec3 origin_set, glm::dvec3 dir_set) 
    {
        ray_origin = origin_set;
        ray_direction = glm::normalize(dir_set);
    }

    //data section
    glm::dvec3 ray_origin;
    glm::dvec3 ray_direction;
    glm::dvec3 intersect_point;
    glm::dvec3 ray_color;

    Triangle* tri_intersect;
    Sphere* sph_intersect;

    double t_val;
};

//ray data structure
std::vector<std::vector<Ray>> all_rays;

const float aspect_ratio = static_cast<float>(WIDTH) / HEIGHT;

unsigned char buffer[HEIGHT][WIDTH][3];

void Arr2Vec(const double double_arr[], glm::dvec3& vec)
{
    vec = glm::dvec3(double_arr[0], double_arr[1], double_arr[2]);
}

/// <summary>
/// Contains Barycentric method of testing to see if a point is contained in a triangle
/// </summary>
/// <param name="tri_check"> -- the triangle to check </param> 
/// <param name="point_to_check" -- the point to check></param>
/// <returns></returns>
bool DoesRayIntersectTriangle(const Triangle* tri_check, glm::dvec3 point_to_check)
{
    //get the three vertices of the triangle
    glm::dvec3 a_vertex, b_vertex, c_vertex;

    Arr2Vec(tri_check->v[0].position, a_vertex);
    Arr2Vec(tri_check->v[1].position, b_vertex);
    Arr2Vec(tri_check->v[2].position, c_vertex);

    //first, ensure that the intersection happens inside the triangle
    glm::dvec3 first_vector = c_vertex - a_vertex;
    glm::dvec3 second_vector = b_vertex - a_vertex;

    //get all dots
    double first_d_point = dot(first_vector, point_to_check);
    double second_d_point = dot(second_vector, point_to_check);
    double first_d_first = dot(first_vector, first_vector);
    double first_d_second = dot(first_vector, second_vector);
    double second_d_second = dot(second_vector, second_vector);

    // barycentric calc
    double area = (first_d_first * second_d_second - first_d_second * first_d_second);
    double a_val = (second_d_second * first_d_point - first_d_second * second_d_point) / area;
    double b_val = (first_d_first * second_d_point - first_d_second * first_d_point) / area;

    // Check if point is in triangle
    return ((a_val + b_val < 1.0) && (a_val >= 0) && (b_val >= 0));
}

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);


/// <summary>
/// Calculates using the barycentric coordinate method
/// </summary>
/// <param name="triangle"></param>
/// <param name="curr_point" - the point the examine></param>
/// <param name="bary_c_val" - output parameter that will contain the barycentric values></param>
void GetBarycentricVal(const Triangle& triangle, const glm::dvec3 curr_point, glm::dvec3& bary_c_val)
{
    glm::dvec3 a_vert, b_vert, c_vert;

    Arr2Vec(triangle.v[0].position, a_vert);
    Arr2Vec(triangle.v[1].position, b_vert);
    Arr2Vec(triangle.v[2].position, c_vert);

    glm::dvec3 first_side = b_vert - a_vert;
    glm::dvec3 second_side = c_vert - a_vert;
    glm::dvec3 vec_to_point = curr_point - a_vert;

    double point_d_first = dot(vec_to_point, first_side);
    double point_d_second = dot(vec_to_point, second_side);
    double first_d_first = dot(first_side, first_side);
    double first_d_sec = dot(first_side, second_side);
    double second_d_second = dot(second_side, second_side);

    //get the dot calculation
    double calc = first_d_first * second_d_second - first_d_sec * first_d_sec;

    //get alpha
    bary_c_val.y = (second_d_second * point_d_first - first_d_sec * point_d_second) / calc;

    //get beta
    bary_c_val.z = (first_d_first * point_d_second - first_d_sec * point_d_first) / calc;

    //get gamma
    bary_c_val.x = 1.0 - bary_c_val.z - bary_c_val.y;
}


/// <summary>
/// Calculates intersection between a ray and a sphere. 
/// if intersection takes place, the ray is updated with the correct
/// intersection values and sphere pointer
/// </summary>
/// <param name="in_sphere"></param>
/// <param name="ray"></param>
void RayIntersectSwphere(const Sphere* in_sphere, Ray& ray) 
{
    for (int i = 0; i < num_spheres; i++) 
    {
        Sphere* this_sphere = &spheres[i];

        if (this_sphere != in_sphere) 
        {
            glm::dvec3 sphere_pos;
            double new_t_val;

            Arr2Vec(this_sphere->position, sphere_pos);

            double b_calc = 2 * (ray.ray_direction.x * (ray.ray_origin.x - sphere_pos.x) + 
                            ray.ray_direction.y * (ray.ray_origin.y - sphere_pos.y) + 
                            ray.ray_direction.z * (ray.ray_origin.z - sphere_pos.z));

            double c_calc = pow((ray.ray_origin.x - sphere_pos.x),2) + 
                            pow((ray.ray_origin.y - sphere_pos.y),2) + 
                            pow((ray.ray_origin.z - sphere_pos.z),2) - 
                            pow(this_sphere->radius,2);

            double b_plus = (-b_calc + sqrt(pow(b_calc, 2) - 4 * c_calc)) / 2.0;
            double b_minus = (-b_calc - sqrt(pow(b_calc, 2) - 4 * c_calc)) / 2.0;

            if (b_plus > b_minus)
            {
                if (b_minus > 0) new_t_val = b_minus;
                else if (b_plus > 0) new_t_val = b_plus;
                else new_t_val = -100.0;
            }
            else
            {
                if (b_plus > 0) new_t_val = b_plus;
                else if (b_minus > 0) new_t_val = b_minus;
                else new_t_val = -100.0;
            }

            if (new_t_val > 0 && (ray.t_val > new_t_val || ray.sph_intersect == nullptr))
            {
                ray.sph_intersect = this_sphere;
                ray.t_val = new_t_val;
                ray.intersect_point = ray.ray_origin + new_t_val * ray.ray_direction;
            }
        }
    }
}


/// <summary>
/// checking to see if the ray intersects with a triangle
//  if so, then update the rays triangle pointer, t value of when intersection occurs
//  and the point on the triangle that intersection takes place
/// </summary>
/// <param name="triangle" a pointer to a triangle that the ray intersected with, if no intersection; nullptr></param>
/// <param name="curr_ray" the ray to be calculated></param>
void RayIntersectTriangle(const Triangle* in_triangle, Ray& curr_ray) 
{
    //loop through all thriangles...
    for (int i = 0; i < num_triangles; i++) 
    {
        Triangle* this_triangle = &triangles[i];

        //...ensuring not to check the trinagle the ray already intersected with
        // if no intersection, [triangle] will be nullptr
        if (this_triangle != in_triangle) 
        {
            //get any arbitrary vertex of a triangle
            glm::dvec3 arb;

            Arr2Vec(this_triangle->v[0].position, arb);

            //ensuring that there is no orthogonallity
            if (dot(this_triangle->tri_normal, curr_ray.ray_direction) != 0) 
            {
                //getting a potential t value by dotting the triangles normal with a vector from an arbitrary point on the triangle
                //and the origin of the ray, dividing it by the angle between the normal and the direction of the ray.
                double new_t_val = dot(this_triangle->tri_normal, arb - curr_ray.ray_origin) / (dot(this_triangle->tri_normal, curr_ray.ray_direction));

                //if the t is in "front" of the ray.
                if (new_t_val > 0) 
                {
                    //if there is no current intersection or if this current t value is closer to the ray's origin
                    if (curr_ray.t_val > new_t_val || curr_ray.tri_intersect == nullptr)
                    {
                        if(DoesRayIntersectTriangle(this_triangle, curr_ray.ray_origin + new_t_val * curr_ray.ray_direction - arb))
                        {
                            curr_ray.tri_intersect = this_triangle;
                            curr_ray.t_val = new_t_val;
                            curr_ray.intersect_point = curr_ray.ray_origin + new_t_val * curr_ray.ray_direction;
                        }
                    }
                }
            }
        }
    }
}


/// <summary>
/// Iterate through all light sources, firing shadow rays from points of intersection on
/// the ray
/// </summary>
/// <param name="curr_ray"></param>
void GetShadowRay(Ray& curr_ray) 
{
    for (int i = 0; i < num_lights; i++) 
    {
        Light* curr_light = &lights[i];

        // shadow ray will only be cast when intersection took place
        if (curr_ray.sph_intersect != nullptr || curr_ray.tri_intersect != nullptr) {

            glm::dvec3 light_pos;
            glm::dvec3 light_reflect;

            Arr2Vec(curr_light->position, light_pos);
            Arr2Vec(curr_light->position, light_reflect);

            //get a vector from the point to the light
            light_reflect = light_reflect - curr_ray.intersect_point;

            //normalize it
            light_reflect = normalize(light_reflect);

            //create a shadow ray starting from the point of intersection heading towards the light
            Ray curr_shadow(curr_ray.intersect_point, light_reflect);

            //calculate intersections of spheres and triangles
            RayIntersectSwphere(curr_ray.sph_intersect, curr_shadow);
            RayIntersectTriangle(curr_ray.tri_intersect, curr_shadow);

            //if there is an intersection, ensure that it is withing bounds
            if (curr_shadow.sph_intersect != nullptr || curr_shadow.tri_intersect != nullptr)
            {
                //calculate the vector bewteen the intersection point
                glm::dvec3 distance = curr_ray.intersect_point - curr_shadow.intersect_point;

                double shadow_intct_pt = dot(distance, distance);

                //calculate the vectore between the intersection and the light
                distance = light_pos - curr_ray.intersect_point;

                double intct_pt_to_lght = dot(distance, distance);

                //if the point that the shadow ray intersects with is further than the light, don't consider it blocked
                if (intct_pt_to_lght < shadow_intct_pt)
                {
                    curr_shadow.sph_intersect = nullptr;
                    curr_shadow.tri_intersect = nullptr;
                }
            }

            //else there is no interction and thus 
            // there must be phong lighting
            if (curr_shadow.sph_intersect == nullptr && curr_shadow.tri_intersect == nullptr)
            {
                glm::dvec3 k_dif, k_spec, norm, reflect, image_plane, light_color, light_vec(curr_shadow.ray_direction);

                double alpha = 0.0f;

                Arr2Vec(curr_light->color, light_color);

                image_plane = -1.0 * curr_ray.ray_direction;

                //sphere
                if (curr_ray.sph_intersect != nullptr) 
                {
                    Arr2Vec(curr_ray.sph_intersect->color_diffuse, k_dif);
                    Arr2Vec(curr_ray.sph_intersect->color_specular, k_spec);

                    norm = curr_ray.sph_intersect->GetNormal(curr_ray.intersect_point);

                    alpha = curr_ray.sph_intersect->shininess;
                }

                //triangle
                else if (curr_ray.tri_intersect != nullptr) 
                {
                    Triangle* triangle = curr_ray.tri_intersect;

                    Vertex vertices[3] = {triangle->v[0], triangle->v[1], triangle->v[2]};

                    glm::dvec3 a_norm, b_norm, c_norm, a_diff, b_diff, c_diff, a_spec, b_spec, c_spec, bary_c_val;;

                    Arr2Vec(vertices[0].normal, a_norm);
                    Arr2Vec(vertices[1].normal, b_norm);
                    Arr2Vec(vertices[2].normal, c_norm);
                    Arr2Vec(vertices[0].color_diffuse, a_diff);
                    Arr2Vec(vertices[1].color_diffuse, b_diff);
                    Arr2Vec(vertices[2].color_diffuse, c_diff);
                    Arr2Vec(vertices[0].color_specular, a_spec);
                    Arr2Vec(vertices[1].color_specular, b_spec);
                    Arr2Vec(vertices[2].color_specular, c_spec);
                        
                    GetBarycentricVal(*triangle, curr_ray.intersect_point, bary_c_val);

                    norm = normalize(a_norm * bary_c_val.x + b_norm * bary_c_val.y + c_norm * bary_c_val.z);

                    k_dif = a_diff * bary_c_val.x + b_diff * bary_c_val.y + c_diff * bary_c_val.z;
                    k_spec = a_spec * bary_c_val.x + b_spec * bary_c_val.y + c_spec * bary_c_val.z;

                    alpha = vertices[0].shininess * bary_c_val.x + vertices[1].shininess * bary_c_val.y + vertices[2].shininess * bary_c_val.z;
                }

                double light_d_norm = dot(light_vec, norm);

                reflect = -1.0 * glm::reflect(light_vec, norm);

                double reflect_d_image = dot(reflect, image_plane);

                //clamping section
                if (reflect_d_image < 0)
                {
                    reflect_d_image = 0;
                }

                if (light_d_norm < 0)
                {
                    light_d_norm = 0;
                }

                //color calculation section
                curr_ray.ray_color.r += light_color.x * (k_dif.x * (light_d_norm)+k_spec.x * pow(reflect_d_image, alpha)) * 255.0;
                curr_ray.ray_color.g += light_color.y * (k_dif.y * (light_d_norm)+k_spec.y * pow(reflect_d_image, alpha)) * 255.0;
                curr_ray.ray_color.b += light_color.z * (k_dif.z * (light_d_norm)+k_spec.z * pow(reflect_d_image, alpha)) * 255.0;
            }
        }
    }
}

//MODIFY THIS FUNCTION
void draw_scene()
{
    glPointSize(2.0);
    glBegin(GL_POINTS);

    for (unsigned int x = 0; x < WIDTH; x++)
    {
        for (unsigned int y = 0; y < HEIGHT; y++)
        {
            //get intersections for all rays and get the shadow rays as well.
            RayIntersectTriangle(nullptr, all_rays[x][y]);
            RayIntersectSwphere(nullptr, all_rays[x][y]);
            GetShadowRay(all_rays[x][y]);

            //change the value of the floats in each of the color channels to change the color of the background
            //when no intersection occurs
            if (all_rays[x][y].sph_intersect == nullptr && all_rays[x][y].tri_intersect == nullptr) 
            {
                plot_pixel(x, y, 1.0f * 255, 1.0f * 255, 1.0f * 255);
            }

            //clamp to 0 or 255 before actually plotting the color
            else 
            {
                glm::dvec3 fin_col;

                fin_col.r = std::min(std::max(all_rays[x][y].ray_color.r + 255.0 * ambient_light[0], 0.0), 255.0);
                fin_col.g = std::min(std::max(all_rays[x][y].ray_color.g + 255.0 * ambient_light[0], 0.0), 255.0);
                fin_col.b = std::min(std::max(all_rays[x][y].ray_color.b + 255.0 * ambient_light[0], 0.0), 255.0);

                plot_pixel(x, y, fin_col.r, fin_col.g, fin_col.b);
            }
        }
    }
    glEnd();
    glFlush();

    printf("Done!\n"); fflush(stdout);
}

void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  buffer[y][x][0] = r;
  buffer[y][x][1] = g;
  buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      //get the three vertices of the triangle
      glm::dvec3 a_vertex, b_vertex, c_vertex;

      Arr2Vec(t.v[0].position, a_vertex);
      Arr2Vec(t.v[1].position, b_vertex);
      Arr2Vec(t.v[2].position, c_vertex);

      //get the trinagle's normal
      t.tri_normal = normalize(cross((b_vertex - a_vertex), (c_vertex - a_vertex)));

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);

  double x_delta, y_delta, x_shift, y_shift;

  //screen calculations
  double tangent = tan(fov * M_PI / 180.0 / 2.0);
  double x_max = tangent * aspect_ratio;
  double y_max = tangent;
  double x_min = -tangent * aspect_ratio;
  double y_min = -tangent;

  //set up increment values
  x_delta = (x_max - x_min) / static_cast<float>(WIDTH);
  y_delta = (y_max - y_min) / static_cast<float>(HEIGHT);

  x_shift = x_min;
  y_shift = y_min;

  //initialize the rays
  all_rays.resize(WIDTH);

  for (int i = 0; i < WIDTH; i++)
  {
      all_rays[i].resize(HEIGHT);
  }

  // fill in all rays with appropriate direction vectors
  for (int i = 0; i < WIDTH; i++) {
      for (int j = 0; j < HEIGHT; j++) {
          all_rays[i][j].Set_Ray(zero_vector, glm::dvec3((i * x_delta) + x_min, (j * y_delta) + y_min, -1.0f));
          y_shift = y_shift + y_delta;
      }
      y_shift = y_min;
      x_shift = x_shift + x_delta;
  }
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 3))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
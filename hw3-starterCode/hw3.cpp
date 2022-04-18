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
#define WIDTH 160
#define HEIGHT 120

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
};

struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};

struct Light
{
  double position[3];
  double color[3];
};

struct SceneObj {
    SceneObj() : obj_num(-1) {};
    std::string obj_type;
    int obj_num;
    double obj_val;
    glm::dvec3 intersect_point;
};

struct Ray {

    Ray(): ray_color(zero_vector) {}; 
    
    Ray(glm::dvec3 in_orig, glm::dvec3 in_dir)
    {
        origin = in_orig;
        direction = glm::normalize(in_dir);
        ray_color = zero_vector;
    }

    void Set_Ray(glm::dvec3 origin_set, glm::dvec3 dir_set) 
    {
        origin = origin_set;
        direction = glm::normalize(dir_set);
    }

    glm::dvec3 origin;
    glm::dvec3 direction;

    SceneObj obj_intersected;
    glm::vec3 ray_color;
};

std::vector<std::vector<Ray>> all_rays;

const float aspect_ratio = (float)WIDTH / HEIGHT;

unsigned char buffer[HEIGHT][WIDTH][3];

double quadraticMinimum(double a, double b, double c) {
    double t0 = (-b + sqrt(b * b - 4 * a * c)) / 2;
    double t1 = (-b - sqrt(b * b - 4 * a * c)) / 2;
    if (t0 > t1) {
        if (t1 > 0) {
            return t1;
        }
        else if (t0 > 0) {
            return t0;
        }
        else {
            return -1;
        }
    }
    else {
        if (t0 > 0) {
            return t0;
        }
        else if (t1 > 0) {
            return t1;
        }
        else {
            return -1;
        }
    }
}

glm::dvec3 toVec3(const double* array) 
{
    return glm::vec3(array[0], array[1], array[2]);
}

//taken from GLM website
double clamp(double x, double minVal, double maxVal) {
    return std::min(std::max(x, minVal), maxVal);
}

//The amount to increment the pixels X and Y values for the rays
double deltaX, deltaY;
int shadowsCount = 0, outerloop = 0, elseloop = 0;

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


glm::dvec3 calcBarycentric(glm::dvec3 point, glm::dvec3 a, glm::dvec3 b, glm::dvec3 c) {
    glm::dvec3 v0 = b - a;
    glm::dvec3 v1 = c - a;
    glm::dvec3 v2 = point - a;

    double d00 = dot(v0, v0);
    double d01 = dot(v0, v1);
    double d11 = dot(v1, v1);
    double d20 = dot(v2, v0);
    double d21 = dot(v2, v1);

    double result = d00 * d11 - d01 * d01;
    glm::dvec3 vec_result;
    vec_result.y = (d11 * d20 - d01 * d21) / result; //alpha
    vec_result.z = (d00 * d21 - d01 * d20) / result; //beta
    vec_result.x = 1.0 - vec_result.z - vec_result.y; //gamma
    return vec_result;
}

void calculateRaySphereIntersection(Ray& ray, int num) {
    for (int i = 0; i < num_spheres; i++) {
        if (i != num) {
            double radius = spheres[i].radius;
            double xc = spheres[i].position[0];
            double yc = spheres[i].position[1];
            double zc = spheres[i].position[2];

            double xd = ray.direction.x;
            double yd = ray.direction.y;
            double zd = ray.direction.z;

            double x0 = ray.origin.x;
            double y0 = ray.origin.y;
            double z0 = ray.origin.z;

            double b = 2 * (xd * (x0 - xc) + yd * (y0 - yc) + zd * (z0 - zc));
            double c = pow((x0 - xc),2) + pow((y0 - yc),2) + pow((z0 - zc),2) - pow(radius,2);
            double result = quadraticMinimum(1, b, c);
            if (result > 0) {
                if (ray.obj_intersected.obj_num == -1 || ray.obj_intersected.obj_val > result) {
                    SceneObj this_obj;
                    this_obj.obj_type = "SPHERE";
                    this_obj.obj_num = i;
                    this_obj.obj_val = result;
                    this_obj.intersect_point = ray.origin + result * ray.direction;
                    ray.obj_intersected = this_obj;
                }
            }
        }
    }
}

void calculateRayTriangleIntersection(Ray& ray, int num) {
    for (int i = 0; i < num_triangles; i++) {
        if (i != num) {
            Triangle triangle = triangles[i];
            glm::dvec3 pointA = glm::dvec3(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
            glm::dvec3 pointB = glm::dvec3(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
            glm::dvec3 pointC = glm::dvec3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

            glm::dvec3 n = cross((pointB - pointA), (pointC - pointA));
            n = normalize(n);
            if (dot(n, ray.direction) != 0) {

                double t = dot(n, pointA - ray.origin) / (dot(n, ray.direction));

                if (t > 0) {
                    if (ray.obj_intersected.obj_num == -1 || ray.obj_intersected.obj_val > t) {
                        // Check if the intersection point is inside the triangle.       
                        glm::dvec3 v0 = pointC - pointA;
                        glm::dvec3 v1 = pointB - pointA;
                        glm::dvec3 v2 = ray.origin + t * ray.direction - pointA;

                        // Compute dot products
                        double dot00 = dot(v0, v0);
                        double dot01 = dot(v0, v1);
                        double dot02 = dot(v0, v2);
                        double dot11 = dot(v1, v1);
                        double dot12 = dot(v1, v2);

                        // Compute barycentric coordinates
                        double denom = (dot00 * dot11 - dot01 * dot01);
                        double u = (dot11 * dot02 - dot01 * dot12) / denom;
                        double v = (dot00 * dot12 - dot01 * dot02) / denom;

                        // Check if point is in triangle
                        if ((u >= 0) && (v >= 0) && (u + v < 1)) {

                            SceneObj this_obj;
                            this_obj.obj_type = "TRIANGLE";
                            this_obj.obj_num = i;
                            this_obj.obj_val = t;
                            this_obj.intersect_point = ray.origin + t * ray.direction;
                            ray.obj_intersected = this_obj;
                        }
                    }
                }
            }
        }
    }
}

/**
 * Fire a shadow ray, and then calculate the color using Phong illumination model if there is no intersection
 * If the shadow ray is obstructed, then return black
 *
 *
 *
 */
void calculateShadowRay(Ray& ray) {
    //Fire a shadow ray first for each light source
    for (int i = 0; i < num_lights; i++) {
        Light light = lights[i];

        //If the ray actually intersected with something, fire the shadow ray
        if (ray.obj_intersected.obj_num != -1) {
            glm::dvec3 lightVec = glm::dvec3(light.position[0], light.position[1], light.position[2]);
            lightVec -= ray.obj_intersected.intersect_point;
            lightVec = normalize(lightVec);
            Ray shadowRay = Ray(ray.obj_intersected.intersect_point, lightVec);

            if (ray.obj_intersected.obj_type == "SPHERE") {
                calculateRaySphereIntersection(shadowRay, ray.obj_intersected.obj_num);
                calculateRayTriangleIntersection(shadowRay, -1);
            }
            else {
                calculateRayTriangleIntersection(shadowRay, ray.obj_intersected.obj_num);
                calculateRaySphereIntersection(shadowRay, -1);
            }

            //Check if shadow ray intersection, if it exists, is behind the light or not
            if (shadowRay.obj_intersected.obj_num != -1) {

                //calculate the vector bewteen the intersection point
                glm::dvec3 distance = ray.obj_intersected.intersect_point - shadowRay.obj_intersected.intersect_point;
                double distanceFromPointToIntersection = dot(distance, distance);

                //calculate the vectore between the intersection and the light
                distance = glm::dvec3(light.position[0], light.position[1], light.position[2]) - ray.obj_intersected.intersect_point;

                double distanceFromPointToLight = dot(distance, distance);

                //if the point that the shadow ray intersects with is further than the light, don't consider it blocked
                if (distanceFromPointToIntersection > distanceFromPointToLight) {
                    shadowRay.obj_intersected.obj_num = -1;
                }
            }

            //if there is no intersection, calculate color using Phong Illumination model with respect to that light
            if (shadowRay.obj_intersected.obj_num == -1) {
                glm::dvec3 kd, ks;
                double alpha = 0.0f; //diffuse, specular, and alpha (shininess)

                glm::dvec3 l, n, r, v, L; //Light vector, normal vector, reflect vector, vector to image plane, Light color
                L = toVec3(light.color);
                v = -ray.direction;
                l = normalize(shadowRay.direction);

                //Calculate lighting for Spheres
                if (ray.obj_intersected.obj_type == "SPHERE") {
                    Sphere s = spheres[ray.obj_intersected.obj_num];

                    //calculate the normal using the vector [from the centery to the intersection] divided by the sphere's radius
                    n = (ray.obj_intersected.intersect_point - glm::dvec3(s.position[0], s.position[1], s.position[2])) / s.radius;

                    kd = glm::dvec3(s.color_diffuse[0], s.color_diffuse[1], s.color_diffuse[2]);
                    ks = glm::dvec3(s.color_specular[0], s.color_specular[1], s.color_specular[2]);
                    alpha = s.shininess;

                }

                //Calculate lighting for triangles
                else if (ray.obj_intersected.obj_type == "TRIANGLE") {
                    Triangle t = triangles[ray.obj_intersected.obj_num];
                    Vertex a = t.v[0], b = t.v[1], c = t.v[2];
                    glm::dvec3 bary = calcBarycentric(ray.obj_intersected.intersect_point, toVec3(a.position), toVec3(b.position), toVec3(c.position));
                    n = normalize(toVec3(a.normal) * bary.x + toVec3(b.normal) * bary.y + toVec3(c.normal) * bary.z);

                    kd = toVec3(a.color_diffuse) * bary.x + toVec3(b.color_diffuse) * bary.y + toVec3(c.color_diffuse) * bary.z;
                    ks = toVec3(a.color_specular) * bary.x + toVec3(b.color_specular) * bary.y + toVec3(c.color_specular) * bary.z;
                    alpha = a.shininess * bary.x + b.shininess * bary.y + c.shininess * bary.z;
                }

                double ln = dot(l, n);
                if (ln < 0) ln = 0;

                r = 2 * (ln)*n - l;
                double rv = dot(r, v);
                if (rv < 0) rv = 0;

                ray.ray_color.r += L.x * (kd.x * (ln)+ks.x * pow(rv, alpha)) * 255;
                ray.ray_color.g += L.y * (kd.y * (ln)+ks.y * pow(rv, alpha)) * 255;
                ray.ray_color.b += L.z * (kd.z * (ln)+ks.z * pow(rv, alpha)) * 255;
            }
        }
    }
}

//MODIFY THIS FUNCTION
void draw_scene()
{
    /*
      //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
      plot_pixel(x, y, x % 256, y % 256, (x+y) % 256);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
    
    */
    //a simple test output
    glPointSize(2.0);
    glBegin(GL_POINTS);

    for (unsigned int x = 0; x < WIDTH; x++)
    {

        for (unsigned int y = 0; y < HEIGHT; y++)
        {
            calculateRayTriangleIntersection(all_rays[x][y], -1);
            calculateRaySphereIntersection(all_rays[x][y], -1);
            calculateShadowRay(all_rays[x][y]);

            //if you can't find the intersection, plot a white color
            if (all_rays[x][y].obj_intersected.obj_num == -1) {
                plot_pixel(x, y, 1.0f * 255, 1.0f * 255, 1.0f * 255);
            }

            //plot the actual color of the ray intersection
            else {
                double r = clamp(all_rays[x][y].ray_color.r + 255 * ambient_light[0], 0, 255);
                double g = clamp(all_rays[x][y].ray_color.g + 255 * ambient_light[1], 0, 255);
                double b = clamp(all_rays[x][y].ray_color.b + 255 * ambient_light[2], 0, 255);

                plot_pixel(x, y, r, g, b);
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

  //Draw Rays
  /*
    rays = new Ray * [WIDTH];
  for (int i = 0; i < WIDTH; i++) {
      rays[i] = new Ray[HEIGHT];
  }
  */

  //initialize the rays
  all_rays.resize(WIDTH);

  for (int i = 0; i < WIDTH; i++)
  {
      all_rays[i].resize(HEIGHT);
  }

  //define the four corners
  double tangentValue = tan(fov * M_PI / 180.0 / 2);
  double x_max = aspect_ratio * tangentValue;
  double x_min = -1 * aspect_ratio * tangentValue;
  double y_max = tangentValue;
  double y_min = -1 * tangentValue;

  /*
    all_rays[0][0] = Ray(0, 0, 0, x_min, y_min, -1);
  all_rays[WIDTH - 1][0] = Ray(0, 0, 0, x_max, y_min, -1);
  all_rays[0][HEIGHT - 1] = Ray(0, 0, 0, x_min, y_max, -1);
  all_rays[WIDTH - 1][HEIGHT - 1] = Ray(0, 0, 0, x_max, y_max, -1);
  */

  all_rays[0][0].Set_Ray(zero_vector, glm::dvec3(x_min, y_min, -1.0f));
  all_rays[WIDTH - 1][0].Set_Ray(zero_vector, glm::dvec3(x_max, y_min, -1.0f));
  all_rays[0][HEIGHT - 1].Set_Ray(zero_vector, glm::dvec3(x_min, y_max, -1.0f));
  all_rays[WIDTH - 1][HEIGHT - 1].Set_Ray(zero_vector, glm::vec3(x_max, y_max, -1.0f));

  //set up increment values
  deltaX = (x_max - x_min) / WIDTH;
  deltaY = (y_max - y_min) / HEIGHT;

  double x_count = x_min;
  double y_count = y_min;

  //create the remaining rays
  for (int i = 0; i < WIDTH; i++) {
      for (int j = 0; j < HEIGHT; j++) {
          all_rays[i][j].Set_Ray(zero_vector, glm::dvec3(x_min + i * deltaX, y_min + j * deltaY, -1));
          y_count += deltaY;
      }
      y_count = y_min;
      x_count += deltaX;
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


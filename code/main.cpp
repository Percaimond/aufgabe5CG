#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <tuple>
#include <vector>

#include "pointlight.h"
#include "scene.h"
#include "sceneobject.h"
#include "util.h"
#include "vec3.h"

const static int SEED = 42;

const static int WIDTH = 600;
const static int HEIGHT = 600;
const static int MAX_DEPTH = 5;

//////////
// TODO 2:
// Compute Phong lighting
//
Vec3d computePhongLighting(
    Vec3d const& view_direction,            //< direction from surface point to camera origin
    Vec3d const& surface_normal,            //< normal vector at surface point
    Vec3d const& light_direction,           //< direction from surface point to light source
    PhongCoefficients const& phong_coeff,   //< phong coefficient k_a,k_d,k_s and n
    Vec3d const& light_color,               //< color of the light source
    double light_intensity)                 //< intensity of the light source
{
    //////////
    // TODO 2:
    // Compute phong lighting.

    Vec3d I = std::get<0>(phong_coeff) * light_intensity;
    Vec3d L = light_direction;

    double NdotL = surface_normal.dot(L);
    if(NdotL > 0) {
        I += std::get<1>(phong_coeff) * light_intensity * NdotL;
        
        Vec3d R = 2. * surface_normal * NdotL - L;
        R = R.normalize();

        double RdotV = R.dot(view_direction);
        if(RdotV > 0) {
            I += std::get<2>(phong_coeff) * light_intensity * std::pow(RdotV, std::get<3>(phong_coeff));
        }
    }
    //I = I.normalize();

    //std::cout << I << std::endl;

    return I; // replace
    // END TODO 2
    /////////////

}

/**
 * @brief Method to check a ray for intersections with any object of the scene.
 * @param ray The ray to trace.
 * @param objects Vector of pointers to all scene objects.
 * @param t_near The intersection distance from the ray origin to the closest point hit.
 * @param hitObject The closest object hit.
 * @return true on hit, false otherwise
 */
bool trace(const Ray& ray,
    const std::vector<std::shared_ptr<SceneObject>>& objects,
    double& t_near, std::shared_ptr<SceneObject>& hitObject)
{
    t_near = std::numeric_limits<double>::max();

    // Check all objects if they got hit by the traced ray. (cf. lecture slide 54)
    // If any object got hit, return the one closest to the camera as 'hitObject'.
    for (auto& o : objects)
    {
        double t = std::numeric_limits<double>::max();

        if (o->intersect(ray, t) && t < t_near)
        {
            hitObject = o;
            t_near = t;
        }
    }

    return (hitObject != nullptr);
}

/**
 * @brief Cast a ray into the scene. If the ray hits at least one object,
 *        the color of the object closest to the camera is returned.
 * @param ray The ray that's being cast.
 * @param objects All scene objects.
 * @return The color of a hit object that is closest to the camera.
 *         Return dark blue if no object was hit.
 */
Vec3d castRay(const Ray& ray, const std::vector<std::shared_ptr<SceneObject>>& objects,
    const std::vector<Pointlight>& lights)
{
    // set the background color as dark blue
    Vec3d hitColor(0, 0, 0.2);

    // early exit if maximum recursive depth is reached - return background color
    if (ray.depth > MAX_DEPTH)
        return hitColor;

    // pointer to the object that was hit by the ray
    std::shared_ptr<SceneObject> hitObject = nullptr;

    // intersection distance from the ray origin to the point hit
    double t = std::numeric_limits<double>::max();

    // Trace the ray. If an object gets hit, calculate the hit point and
    // retrieve the surface color 'hitColor' from the 'hitObject' object that was hit
       if (trace(ray, objects, t, hitObject))
    {
        hitColor = Vec3d();
		// Intersection point with the hit object
        Vec3d p_hit = ray.origin + ray.dir * t;

        for (unsigned i=0; i < lights.size(); i++) {
            Ray shadow_ray;
            shadow_ray.origin = p_hit;
            shadow_ray.dir = lights[i].getPosition();
            double nearest;
            std::shared_ptr<SceneObject> shadow_hit = nullptr;

            if(trace(shadow_ray,objects,nearest,shadow_hit)){
                Vec3d k_a = std::get<0>((*hitObject).getPhongCoefficients(p_hit));
                double light_intensity = 0.01;
                hitColor = k_a * light_intensity;
                //apply ambient term only
                //beleuchten Sie den Schnittpunkt direkt nur mit dem ambienten Term des Phong-Belechtungsmodells
            };

            /*
            Vec3d const& view_direction,			//< direction from surface point to camera origin
            Vec3d const& surface_normal,			//< normal vector at surface point
            Vec3d const& light_direction,			//< direction from surface point to light source
            PhongCoefficients const& phong_coeff,	//< phong coefficient k_a,k_d,k_s and n
            Vec3d const& light_color,				//< color of the light source
            double light_intensity)					//< intensity of the light source
            */

            const Vec3d cameraPos(-5, -1, 0.);
            Vec3d view_direction = Vec3d(cameraPos[0]-p_hit[0], cameraPos[1]-p_hit[1], cameraPos[2]-p_hit[2]);
            Vec3d surface_normal = (*hitObject).getSurfaceNormal(p_hit);
            Vec3d light_pos = lights[i].getPosition();
            Vec3d light_direction = Vec3d(light_pos[0]-p_hit[0], light_pos[1]-p_hit[1], light_pos[2]-p_hit[2]);
            PhongCoefficients phong_coeff = (*hitObject).getPhongCoefficients(p_hit);
            Vec3d light_color = lights[i].getColor();
            double light_intensity = 0.09;
            
            view_direction = (view_direction).normalize();
            surface_normal = (surface_normal).normalize();
            light_direction = (light_direction).normalize();

            hitColor = computePhongLighting(view_direction,surface_normal,light_direction,phong_coeff,light_color,light_intensity);
        }
        //////////
        // TODO 3:
        // Compute local lighting. The result is added to "hitColor".
        //
        // For each light source (given by funtion parameter lights)
        //
        //   a) Cast a shadow ray from the hitpoint to the light-source (use the trace function)
        //
        //   b) If not in shadow, compute local lighting using the function "computePhongLighting"
        //      Else apply ambient term only
        //
        //      For a more realistic image, use inverse square attentuation for the light intensity.
        //

        // END TODO 3
        /////////////

        //////////
        // TODO 4:
        // Compute reflection.
        //
        // Build a reflection for the hitpoint of the hit object. 
        // Use this ray to make a recursive call to "castRay" and add the result of the call to "hitColor".
        //

        // END TODO 4
        /////////////
    }

    return hitColor;
}

/**
 * @brief The rendering method, loop over all pixels in the framebuffer, shooting
 *        a ray through each pixel with the origing being the camera position.
 * @param viewport Size of the framebuffer.
 * @param objects Vector of pointers to all objects contained in the scene.
 */
void render(const Vec3i viewport, const std::vector<std::shared_ptr<SceneObject>>& objects,
    const std::vector<Pointlight>& lights)
{
    std::vector<Vec3d> framebuffer(static_cast<size_t>(viewport[0]) * viewport[1]);

    // camera position in world coordinates
    const Vec3d cameraPos(0., 0., 0.);

    // view plane parameters
    const double l = -1.;   // left
    const double r = +1.;   // right
    const double b = -1.;   // bottom
    const double t = +1.;   // top
    const double d = +2.;   // distance to camera

    // Cast a ray from 'cameraPos' through the center(!) of each pixel on the viewplane.
    // Use the view plane parametrization given above (l,r,b,t,d).
    #pragma omp parallel for
    for (int j = 0; j < viewport[1]; ++j)
    {
        for (int i = 0; i < viewport[0]; ++i)
        {
            double u = l + (r - l) * (i + 0.5) / viewport[0];
            double v = t + (b - t) * (j + 0.5) / viewport[1];

            Ray ray;
            ray.origin = cameraPos;
            ray.dir = Vec3d(u, v, -d) - cameraPos;
            ray.dir = ray.dir.normalize();
            framebuffer.at(i + j * static_cast<size_t>(viewport[0])) = castRay(ray, objects, lights);
        }
    }

    // save the framebuffer an a PPM image
    saveAsPPM("./result.ppm", viewport, framebuffer);
}

/**
 * @brief main routine.
 *        Generates the scene and invokes the rendering.
 * @return
 */
int main()
{
    // Generate the scene objects
    const auto objects = create_scene_objects();

    // Let there be light
    const auto lights = create_scene_lights();

    // Start rendering
    const Vec3i viewport(WIDTH, HEIGHT, 0);
    render(viewport, objects, lights);

    return 0;
}

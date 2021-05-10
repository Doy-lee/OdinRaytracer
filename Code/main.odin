package main

import "core:fmt"
import "core:mem"
import "core:os"
import "core:math"
import "core:math/linalg"

Dqn_V2i :: [2]int;
Dqn_V3f :: [3]f32;
Dqn_V2f :: [2]f32;

BMPHeader :: struct #packed
{
    file_type: u16, // BM in hex for Windows NT etc
    file_size: u32, // Size of the BMP file in bytes
    reserved01: u16, // Application reserved, set to 0 if manually created
    reserved02: u16, // Application reserved, set to 0 if manually created
    offset_to_pixels: u32, // Byte offset to the pixel data in the image

    // DIB Header
    header_size: u32, // Must be 40 for BITMAPINFOHEADER
    width: i32,
    height: i32,
    color_planes: u16,
    bits_per_pixel: u16,
    compression_method: u32, // Specific magic bytes, or 0 for no compression.
    image_size: u32, // The size of the raw bitmap data. Can be 0 if compression is BI_RGB (0)
    horizontal_resolution: u32, // Horizontal pixels per meter
    vertical_resolution: u32, // Veritcal pixels per meter
    colors_in_color_palette_count: u32, // Or set to 0 to default to 2^n
    important_colors_count: u32, // Or set to 0 when every color is important; generally ignored
};

ImageRGBA :: struct
{
    width, height: int,
    pixels: []u32,
};

LightType :: enum{Point, Directional, Ambient};
Light :: struct
{
    type : LightType,
    intensity : f32,
    xyz : Dqn_V3f,
};

Material :: struct
{
    color: Dqn_V3f,
    specular_exponent : f32,
};

Plane :: struct
{
    n: Dqn_V3f,
    d: f32,
    material_index: u32,
};

Sphere :: struct
{
    p: Dqn_V3f,
    r: f32,
    material_index: u32,
};

World :: struct
{
    materials: []Material,
    planes: []Plane,
    spheres: []Sphere,
    lights: []Light,
};

BMPHeader_ImageSize :: proc(header: ^BMPHeader) -> Dqn_V2i
{
    result := Dqn_V2i{cast(int)header.width, cast(int)abs(header.height)};
    return result;
}

BMPHeader_Init :: proc(width, height: i32, bits_per_pixel: u16) -> BMPHeader
{
    pixels_size_in_bytes := cast(int)(width * height * cast(i32)bits_per_pixel / 8);

    result : BMPHeader;
    result.file_type             = 'M' << 8 | 'B' << 0;
    result.file_size             = cast(u32)(size_of(result) + pixels_size_in_bytes);
    result.offset_to_pixels      = size_of(result);
    result.header_size           = 40; // Fixed value for BITMAPINFOHEADER which we support
    result.width                 = width;
    result.height                = height;
    result.color_planes          = 1;
    result.bits_per_pixel        = bits_per_pixel;
    result.horizontal_resolution = 1;
    result.vertical_resolution   = 1;
    return result;
}

ImageRGBA_Write :: proc(image : ImageRGBA, filename: string) -> b32
{
    header := BMPHeader_Init(width = cast(i32)image.width, height = cast(i32)image.height, bits_per_pixel = 32);
    file_handle, error := os.open(filename, os.O_WRONLY | os.O_CREATE);
    result := cast(b32)(error == os.ERROR_NONE);
    if result {
        os.write(file_handle, mem.ptr_to_bytes(&header));
        os.write(file_handle, transmute([]u8)image.pixels);
        os.close(file_handle);
    }
    return result;
}

ImageRGBA_Allocate :: proc(width, height: int) -> ImageRGBA
{
    pixels_size_in_bytes := cast(int)(height * width * 4 /*bpp*/);

    result : ImageRGBA;
    result.width = width;
    result.height = height;
    result.pixels = make([]u32, pixels_size_in_bytes);
    return result;
}

main :: proc()
{
    materials : []Material =
    {
        {color = Dqn_V3f{0, 0, 0}, specular_exponent = 1.0},
        {color = Dqn_V3f{1, 0, 0}, specular_exponent = 1.0},
        {color = Dqn_V3f{0, 0, 1}, specular_exponent = 1.0},
        {color = Dqn_V3f{0, 1, 0}, specular_exponent = 1.0},
    };

    // Equations
    //
    // Ray
    //     Starting from a point 'p1' combined with a direction (p2 - p1) and
    //     't' to move along the vector's direction
    //
    //     p1 + (p2 - p1)t
    //
    // Sphere
    //     A point is considered on the sphere IFF the distance between the
    //     centre of the sphere and the point 'p1' is equal to the radius of the
    //     sphere. All the points that lie on the sphere constitute the equation
    //     of a sphere.
    //
    //     distance(p1, sphere_centre)                       = radius
    //     sqrt(dot(p1 - sphere_centre, p1 - sphere_centre)) = radius
    //     dot(p1 - sphere_centre, p1 - sphere_centre)       = radius^2
    //     (p1 - sphere_centre)^2 - radius^2                 = 0
    //
    //     And, if 'p1 - sphere_centre' was a singular point 'p' then
    //
    //     dot(p, p)             = radius^2
    //     p.x^2 + p.y^2 + p.z^2 = radius^2
    //     p^2                   = radius^2
    //     p^2 - radius^2        = 0
    //
    // Ray intersects sphere surface
    //
    //     Ray               = p1 + (p2 - p1)t
    //     Point on a sphere = p^2 - radius^2 = 0
    //
    //     p^2 = (p1 + (p2 - p1)t) * (p1 + (p2 - p1)t)
    //         = p1^2 + p1(p2 - p1)t + p1(p2 - p1)t + ((p2 - p1)t)^2
    //         = p1^2 + 2*p1(p2 - p1)t + (p2 - p1)^2*t^2
    //
    //     Then 'p^2 - radius^2' becomes
    //
    //     p1^2 + 2*p1(p2 - p1)t + (p2 - p1)^2*t^2 - radius^2 = 0
    //
    //     Re-arrange the equation, focusing on t
    //
    //     (p2-p1)^2*t^2 + 2*p1*(p2-p1)*t + p1^2 - radius^2 = 0
    //              ax^2 +             bx + c               = 0
    //
    //     'p2-p1' is actually the direction of the ray, so let 'p1-p2' be 'dir'
    //
    //     dir^2*t^2 + 2*p1*dir*t + p1^2 - radius^2 = 0
    //
    //     Or in other words, this is a quadratic polynomial and can be solved
    //     with the quadratic equation.
    //
    //     t = (-b ± sqrt(b^2 - 4ac)) / 2a
    //
    //     Note that any squared components is equivalent to the dot product of
    //     it self and any vector multiply is a dot product
    //     (since 'p^2 - radius^2 = 0' is equivalent to 'dot(p, p) - radius^2')
    //
    //     dot(dir)*t^2 + 2*dot(p1, dir)*t + dot(p1, p1) - radius^2 = 0

    planes : []Plane =
    {
        {n = Dqn_V3f{0, 0, 1}, d = 0, material_index = 1},
    };

    spheres : []Sphere =
    {
        {p = Dqn_V3f{0, -1, 3}, r = 1.0, material_index = 1},
        {p = Dqn_V3f{2, 0, 4}, r = 1.0, material_index = 2},
        {p = Dqn_V3f{-2, 0, 4}, r = 1.0, material_index = 3},
    };

    lights : []Light =
    {
        {type = LightType.Ambient, intensity = 0.2},
        {type = LightType.Point, intensity = 0.6, xyz = Dqn_V3f{2, 1, 0}},
        {type = LightType.Directional, intensity = 0.2, xyz = Dqn_V3f{1, 4, 4}},
    };

    world : World;
    world.materials = materials[:];
    world.planes = planes[:];
    world.spheres = spheres[:];
    world.lights = lights[:];

    camera_p := Dqn_V3f{0, 0, 0};

    film_dist : f32 = 1.0;
    film_size := Dqn_V2f{1.0, 1.0};
    film_p : [2]f32 = {(-0.5 * film_size.x), (-0.5 * film_size.y)};

    image : ImageRGBA = ImageRGBA_Allocate(720, 720);
    for y := 0; y < image.height; y +=1
    {
        for x := 0; x < image.width; x +=1
        {
            pixel_to_film_p := Dqn_V3f{film_p.x + (film_size.x / cast(f32)image.width) * cast(f32)x,
                                       film_p.y + (film_size.y / cast(f32)image.height) * cast(f32)y,
                                       film_dist};

            ray_origin    : Dqn_V3f = camera_p;
            ray_direction : Dqn_V3f = linalg.vector_normalize(pixel_to_film_p - camera_p);

            best_t: f32 = math.F32_MAX;
            best_sphere: Sphere;
            for _, i in world.spheres
            {
                sphere := &world.spheres[i];

                // See: Ray Intersects Sphere Surface comment for more info
                //
                // Equation for intersection
                // dot(dir)*t^2 + 2*dot(p1, dir)*t + dot(p1, p1) - radius^2 = 0
                //
                // Is of the form
                // ax^2 + bx + c = 0 (i.e. the quadratic formula)
                //
                // (a)x^2 = dot(dir)
                // (b)x   = 2 * dot(p1, dir)
                // c      = dot(p1, p1) - radius^2
                //
                // Solved with the quadratic equation
                // x = (-b ± sqrt(b^2 - 4ac)) / 2a
                //
                // Note before that the 'p' (or 'p1') in this case is actually
                // represented as 'p1 - sphere_centre' in the sphere equation
                // 'p^2 - r^2 = 0' (as explained in the equation of a sphere
                // above) then we need to ensure p1 is calculated as such.

                sphere_relative_origin : Dqn_V3f = ray_origin - sphere.p;
                a                : f32 = linalg.vector_dot(ray_direction, ray_direction);
                b                : f32 = 2 * linalg.vector_dot(sphere_relative_origin, ray_direction);
                c                : f32 = linalg.vector_dot(sphere_relative_origin, sphere_relative_origin) - (sphere.r * sphere.r);
                discriminant     : f32 = (b * b) - (4 * a * c);

                // Ray intersection. Note if the discriminant was 0, the
                // quadratic equation simplifies down to one solution, i.e.
                // t_plus and t_minus give the same result.
                //
                // If the discriminant is < 0, then we have an square root of
                // a negative number which is an imaginary number. In this case,
                // we don't have a valid solution so the intersection is
                // ignored.
                if discriminant > 0
                {
                    t_plus  : f32 = (-b + math.sqrt(discriminant)) / 2 * a;
                    t_minus : f32 = (-b - math.sqrt(discriminant)) / 2 * a;

                    if (t_plus > film_dist && t_plus < best_t)
                    {
                        best_t = t_plus;
                        best_sphere = sphere^;
                    }

                    if (t_minus > film_dist && t_minus < best_t)
                    {
                        best_t = t_minus;
                        best_sphere = sphere^;
                    }
                }
            }

            material := &materials[best_sphere.material_index];
            light_coeff : f32;
            for _, light_index in world.lights
            {
                light : ^Light = &world.lights[light_index];
                switch light.type
                {
                    case .Ambient:
                        light_coeff += light.intensity;

                    case .Directional: fallthrough;
                    case .Point:

                        // Light is represented as the intensity of the color
                        // (i.e. a coefficient on the color values on point of
                        // the surface). The intensity is calculated as the
                        // projection of the vector from the surface to the
                        // light against the normal of the surface, otherwise
                        // known as the dot product.
                        //
                        // Note that when the light vector is exactly parallel
                        // to the surface, the intensity is calculated as 1,
                        // (i.e. the dot product of 2 normalized vectors that
                        // are parallel is known as dot(a, b) == 1) hence the
                        // color reflected is at it's most intense.
                        //
                        // Similarly, when (dot(a, b) == 0) the projection of
                        // a onto b, means the light hitting the surface at an
                        // angle of 90 degrees and is modelled as having no
                        // impact on the intensity of the color reflected.
                        //
                        surface_p           : Dqn_V3f = ray_origin + (ray_direction * best_t);
                        light_vector        : Dqn_V3f = light.type == .Directional ? light.xyz : light.xyz - surface_p;
                        surface_n           : Dqn_V3f = linalg.vector_normalize(surface_p - best_sphere.p);
                        surface_n_dot_light : f32     = linalg.vector_dot(surface_n, light_vector);
                        if surface_n_dot_light > 0
                        {
                            light_vector_length : f32 = linalg.vector_length(light_vector);
                            light_coeff += light.intensity * surface_n_dot_light / light_vector_length;
                        }

                        /*
                            Calculating a reflection vector given a normal and an incident ray.

                                 N
                             L   ^   R
                             ^   |   ^
                              \  |  /
                               \a|b/
                            ____\|/____
                                 x

                            L Light vector (i.e. incident ray)
                            N Normal to the surface
                            R Reflection vector [reflection of the light vector]
                            a Angle between L and N
                            b Angle between N and R
                            x Intersect point of the light vector (indicent ray) against the surface

                            The reflection vector is defined where angle a == b.
                            Since angle (a == b) then the vector LN is the same length as NR by
                            symmetry.

                                 N
                             L   ^   R
                             ^---n   ^
                              \  |  /
                               \a|b/
                            ____\|/____
                                 x

                            To get the length Ln notated by the '---' in the diagram, it can be
                            calculated as length(L - nx), but we don't know what 'nx' is exactly. We
                            have a normal vector (N) that is of unit length that we can use the dot
                            product on to figure out the length of 'nx' and scale the normal vector
                            via the length.

                            Hence,
                             - Project L onto N via the dot product dot(L, N), this gives us the
                               length of the vector nx
                             - Multiply the length of the vector nx agains the normal vector (N)

                            i.e. nx = dot(L, N) * N
                                    = length(nx) * N

                            then

                                 Ln = L - nx

                            Similarly,

                                 Rn = R - nx

                            But note that Rn is the same as -Ln due to symmetry (a == b), hence

                                 Rn              = -Ln
                                 R - nx          = -(L - nx)
                                 R - dot(L, N)*N = -L + dot(L, N)*N
                                 R               = 2*dot(L, N)*N - L
                         */

                        reflected_light_v3 : Dqn_V3f = (2 * surface_n_dot_light * surface_n) - light_vector;
                        surface_to_ray_origin_v3 : Dqn_V3f = ray_origin - surface_p;

                        reflected_light_v3_dot_surface_to_ray_origin_v3 : f32 = linalg.dot(reflected_light_v3, surface_to_ray_origin_v3);
                        if reflected_light_v3_dot_surface_to_ray_origin_v3 > 0 {
                            reflected_light_v3_length       : f32 = linalg.vector_length(reflected_light_v3);
                            surface_to_ray_origin_v3_length : f32 = linalg.vector_length(surface_to_ray_origin_v3);

                            angle_between_reflected_light_and_ray_origin : f32 = reflected_light_v3_dot_surface_to_ray_origin_v3 / (reflected_light_v3_length * surface_to_ray_origin_v3_length);
                            light_coeff += light.intensity * math.pow(angle_between_reflected_light_and_ray_origin, material.specular_exponent);
                        }
                }
            }

            r := cast(u8)(cast(f32)material.color[0] * 255.0 * light_coeff);
            g := cast(u8)(cast(f32)material.color[1] * 255.0 * light_coeff);
            b := cast(u8)(cast(f32)material.color[2] * 255.0 * light_coeff);
            image.pixels[(y * image.width) + x] = cast(u32)0xff << 24 | cast(u32)r << 16 | cast(u32)g << 8 | cast(u32)b << 0;
        }
    }

    FILE :: "test.bmp";
    if ImageRGBA_Write(image, FILE)
    {
        fmt.printf(`[INFO] Image written to "{}"`, FILE);
    }
    else
    {
        fmt.printf(`[ERROR] Failed to write image to {}`, FILE);
    }
}

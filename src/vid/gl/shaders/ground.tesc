#version 450 core

layout (vertices = 3) out;

uniform vec2  u_res;
uniform mat4  u_mvp;
uniform float u_lod;

#if 1
// control:
//     layout(vertices = 4) out;
// 
//     uniform vec2 screen_size;
//     uniform mat4 mvp;
//     uniform float lod_factor;
    
bool offscreen(vec4 vertex)
{
  if(vertex.z < -0.5) return true;
  return any(lessThan(vertex.xy, vec2(-1.7))) ||
    any(greaterThan(vertex.xy, vec2(1.7)));
  // return any(
  //     lessThan(vertex.xy, vec2(-1.7)) ||
  //     greaterThan(vertex.xy, vec2(1.7))
  //     );
}

vec4 project(vec4 vertex)
{
  vec4 result = u_mvp * vertex;
  result /= result.w;
  return result;
}

vec2 screen_space(vec4 vertex)
{
  return (clamp(vertex.xy, -1.3, 1.3)+1) * (u_res*0.5);
}

float level(vec2 v0, vec2 v1)
{
  // TODO: also clamp to the world space distance, no need to go smaller than texture size
  // this is with the 16x16 characters:
  // const float k_terrain_pixel_scale = 1.0f/(16.0*1024.0/(3.0*0.3048*2048.0));
  // const float k_terrain_pixel_scale = 16.0f/(16.0*1024.0/(3.0*0.3048*2048.0));
  // return clamp(distance(v0, v1)/u_lod, 1, 64);
  return clamp(distance(v0, v1)/16.0f, 1, 64);
}

void main()
{
  if(gl_InvocationID == 0)
  {
    vec4 v0 = project(gl_in[0].gl_Position);
    vec4 v1 = project(gl_in[1].gl_Position);
    vec4 v2 = project(gl_in[2].gl_Position);
    // vec4 v3 = project(gl_in[3].gl_Position);

    // if(all(bvec4(offscreen(v0), offscreen(v1), offscreen(v2), offscreen(v3))))
    if(all(bvec3(offscreen(v0), offscreen(v1), offscreen(v2))))
    {
      gl_TessLevelInner[0] = 0;
      // gl_TessLevelInner[1] = 0;
      gl_TessLevelOuter[0] = 0;
      gl_TessLevelOuter[1] = 0;
      gl_TessLevelOuter[2] = 0;
      // gl_TessLevelOuter[3] = 0;
    }
    else
    {
      vec2 ss0 = screen_space(v0);
      vec2 ss1 = screen_space(v1);
      vec2 ss2 = screen_space(v2);
      // vec2 ss3 = screen_space(v3);

      float e0 = level(ss1, ss2);
      float e1 = level(ss0, ss2);
      float e2 = level(ss1, ss0);
      // float e0 = level(ss1, ss2);
      // float e1 = level(ss0, ss1);
      // float e2 = level(ss3, ss0);
      // float e3 = level(ss2, ss3);

      gl_TessLevelInner[0] = mix(e1, e2, 0.5);
      // gl_TessLevelInner[1] = mix(e0, e3, 0.5);
      gl_TessLevelOuter[0] = e0;
      gl_TessLevelOuter[1] = e1;
      gl_TessLevelOuter[2] = e2;
      // gl_TessLevelOuter[3] = e3;
    }
  }
  gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
}
#endif

#if 0
void main(void)
{
  if (gl_InvocationID == 0)
  {
    // determine three edge tessellation levels by distance in texture space.
    // if the distance is so small that we couldn't resolve the finest texture
    // level, we don't care.
    // also, if we're too far from camera (u_pos_ws) we don't care either
    vec3 a = gl_in[0].gl_Position.xyz;
    vec3 b = gl_in[1].gl_Position.xyz;
    vec3 c = gl_in[2].gl_Position.xyz;
    // this is with the 16x16 characters:
    // const float k_terrain_pixel_scale = 1.0f/(16.0*1024.0/(3.0*0.3048*2048.0));
    const float k_terrain_pixel_scale = 16.0f/(16.0*1024.0/(3.0*0.3048*2048.0));
    float e0 = length(b-c);
    float e1 = length(c-a);
    float e2 = length(a-b);
    float m0 = clamp(1.0f + 2*length(0.5f*(b+c)) / 40.0f, 1.0f, 320.0f);
    float m1 = clamp(1.0f + 2*length(0.5f*(a+c)) / 40.0f, 1.0f, 320.0f);
    float m2 = clamp(1.0f + 2*length(0.5f*(a+b)) / 40.0f, 1.0f, 320.0f);
    gl_TessLevelOuter[0] = clamp(e0 / (m0 * k_terrain_pixel_scale), 1, 8);
    gl_TessLevelOuter[1] = clamp(e1 / (m1 * k_terrain_pixel_scale), 1, 8);
    gl_TessLevelOuter[2] = clamp(e2 / (m2 * k_terrain_pixel_scale), 1, 8);
    gl_TessLevelInner[0] = max(
        gl_TessLevelOuter[0],
        max(gl_TessLevelOuter[1], gl_TessLevelOuter[2]));
  }
  gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
}
#endif

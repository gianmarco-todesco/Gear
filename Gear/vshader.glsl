#version 400

uniform mat4 mvp_matrix;
uniform mat4 mv_matrix;
uniform mat3 normal_matrix;

in vec4 a_position;
in vec3 a_normal;
attribute vec2 a_texcoord;

varying vec3 Position;
varying vec3 Normal;

out Data 
{
	vec3 Position;
	vec3 Normal;
	// vec2 TexCoord;
} data;


//! [0]
void main()
{

    Normal = normalize( normal_matrix * a_normal );
	Position = vec3(mv_matrix * a_position);
	// v_texcoord = data.Normal.xy;
	//

    // Calculate vertex position in screen space
    gl_Position = mvp_matrix * a_position;

    // Pass texture coordinate to fragment shader
    // Value will be automatically interpolated to fragments inside polygon faces
    // v_texcoord = a_texcoord;
}
//! [0]

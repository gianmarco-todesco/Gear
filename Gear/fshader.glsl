#version 400


in Data
{
	vec3 Position;
	vec3 Normal;
	vec2 TexCoord;
} data;

uniform sampler2D texture;

varying vec2 v_texcoord;

varying vec3 Position;
varying vec3 Normal;

//! [0]
void main()
{
    // Set fragment color from texture
	

	vec3 n = normalize( Normal );
	vec3 s = normalize( vec3(0.0,2.0,-1.0) - Position );
	vec3 v = normalize( -Position );
	vec3 r = reflect( -s, n );
 
	vec3 ambient = vec3(0.1,0.1,0.1);
 
	float sDotN = max( dot( s, n ), 0.0 );
	vec3 diffuse = vec3(0.3,0.8,0.01) * sDotN; 
 
	vec3 spec = vec3(0.1,0.1,0.1) * pow( max( dot(r,v) , 0.0 ), 120.0 ); 


    gl_FragColor = vec4(ambient + diffuse + spec, 1.0); //  texture2D(texture, v_texcoord);
}
//! [0]


uniform float omega_w;
uniform float invomega;
uniform vec3 c;
uniform sampler2D streamed; // single channel
uniform sampler2D macros;

void main () {
  vec4 u = texture2D(macros, gl_FragCoord);
  vec3 v = u.xyz;

  float u2 =  1 - 1.5 * dot(v, v);
  float dir = dot(c, v);
  float output = omega_w * u.w * (u2 + 3 * dir + 4.5 * dir * dir) +
                 invomega * texture2D(streamed, gl_FragCoord).r;
  gl_FragColor = vec4(output);
}